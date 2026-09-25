"""Axo-dendritic vs axo-axonic input, split with navis (Philipp Schlegel's
implementation of synaptic flow centrality). Python twin of axoaxo.R.

    uv run axoaxo_navis.py --pre AVLP749m --post PAM02 PAM01 --name test
    uv run axoaxo_navis.py --set-token          # new neuprint token -> .env

or from Python:

    from axoaxo_navis import axoaxo
    axoaxo(pre=["AVLP749m"], post=["PAM02", "PAM01"], name="test")

pre / post: bodyids and/or type names (types resolve via neuprint).
Each postsynaptic body is split with navis.split_axon_dendrite(metric=
"synapse_flow_centrality", split="prepost"), which labels nodes axon, dendrite,
linker (the high-flow bridge between them) or cellbodyfiber. Each synapse from
`pre` takes the compartment of its nearest skeleton node.

Writes to scenes_navis/<name>/: synapses.csv, summary.csv, scene.json,
scene_url.txt and <name>.webloc (double-click to open the scene).
"""
import argparse
import getpass
import json
import os
import re
import subprocess
import urllib.parse
import warnings
from pathlib import Path

import navis
import navis.interfaces.neuprint as nvn
import neuprint as neu
import pandas as pd
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

AXO_DIR = Path(__file__).resolve().parent
CACHE = AXO_DIR / "cache_navis"
DATASET = "male-cns:v1.0"
# bright on neuroglancer's black background; same colours as axoaxo.R
DOMAIN_COLS = {"axo-dendritic": "#40C4FF", "axo-axonic": "#FF5252", "axo-linker": "#FFD740"}
# one colour per neuron type: warm for post types, cool for pre types
POST_PAL = ["#ff0000", "#ff9100", "#ff4081", "#ffd600", "#d500f9", "#ff6e40"]
PRE_PAL = ["#00e5ff", "#76ff03", "#2979ff", "#1de9b6", "#b388ff", "#c6ff00"]
DOMAIN_OF = {"dendrite": "axo-dendritic", "axon": "axo-axonic",
             "linker": "axo-linker", "cellbodyfiber": "axo-linker"}

CACHE.mkdir(exist_ok=True)


def token():
    """neuprint token: environment, then .env in this folder, then ~/.Renviron (R users)."""
    if os.environ.get("NEUPRINT_APPLICATION_CREDENTIALS"):
        return os.environ["NEUPRINT_APPLICATION_CREDENTIALS"]
    for f, key in ((AXO_DIR / ".env", "NEUPRINT_APPLICATION_CREDENTIALS"),
                   (Path("~/.Renviron").expanduser(), "NEUPRINT_TOKEN")):
        if f.exists():
            m = re.search(rf'^{key}="?([^"\n]+)', f.read_text(), re.M)
            if m:
                return m.group(1)
    raise RuntimeError("no neuprint token: run `uv run axoaxo_navis.py --set-token`")


def set_token(new=None):
    """Save a neuprint token to .env, replacing any old one; asks for it if not given."""
    new = (new or getpass.getpass("neuprint token (input hidden): ")).strip()
    if not new:
        raise SystemExit("no token entered, .env unchanged")
    env = AXO_DIR / ".env"
    lines = env.read_text().splitlines() if env.exists() else []
    lines = [l for l in lines if not l.startswith("NEUPRINT_APPLICATION_CREDENTIALS=")]
    env.write_text("\n".join(lines + [f"NEUPRINT_APPLICATION_CREDENTIALS={new}"]) + "\n")
    print(f"saved to {env}")
    if os.environ.get("NEUPRINT_APPLICATION_CREDENTIALS"):
        print("note: NEUPRINT_APPLICATION_CREDENTIALS is also set in your shell and takes "
              "precedence; unset it to use the new token")


def client():
    os.environ["NEUPRINT_APPLICATION_CREDENTIALS"] = token()
    for srv in (os.environ.get("NEUPRINT_HOST", "https://neuprint.janelia.org"),
                "https://neuprint-cns.janelia.org"):
        try:
            return neu.Client(srv, dataset=DATASET)
        except Exception:
            pass
    raise RuntimeError("neuprint login failed on both hosts")


def resolve_ids(x, c):
    """Bodyids pass through; anything else is a type name."""
    x = [str(v) for v in x]
    types = [v for v in x if not v.isdigit()]
    hits = neu.fetch_neurons(neu.NeuronCriteria(type=types), client=c)[0] if types else pd.DataFrame(
        columns=["bodyId", "type"])
    miss = set(types) - set(hits.type)
    if miss:
        raise ValueError(f"unknown type(s): {', '.join(sorted(miss))}")
    ids = [int(v) for v in x if v.isdigit()] + hits.bodyId.tolist()
    return list(dict.fromkeys(int(i) for i in ids))


def fetch_skeleton(bid, c):
    # pickle is fine here: the cache only ever holds files this function wrote
    f = CACHE / f"{bid}.pkl"
    if f.exists():
        return pd.read_pickle(f)
    n = nvn.fetch_skeletons(bid, with_synapses=True, heal=True, client=c)[0]
    pd.to_pickle(n, f)
    return n


def split(n, min_sites):
    """navis split -> (node_id -> compartment, segregation index), or None."""
    cn = n.connectors
    if (cn.type == "post").sum() < min_sites or (cn.type == "pre").sum() < min_sites:
        return None
    kw = dict(metric="synapse_flow_centrality", split="prepost", reroot_soma=True)
    lab = navis.split_axon_dendrite(n, label_only=True, **kw)
    seg = navis.segregation_index(navis.split_axon_dendrite(n, **kw))
    return lab.nodes.set_index("node_id").compartment, seg


def mcns_scene(title):
    """Same base scene as neuroglancer_scene.R::make_mcns_scene, minus the
    single malecns layer: segmentation layers are added per type by the caller."""
    shell = lambda url, name, tab: dict(
        type="segmentation",
        source=dict(url=url, subsources=dict(default=True, properties=True, mesh=True),
                    enableDefaultSubsources=False),
        pick=False, tab=tab, selectedAlpha=0, saturation=0, meshSilhouetteRendering=7,
        segments=["1"], colorSeed=1336242844, segmentDefaultColor="#ffffff", name=name)
    roi = "precomputed://gs://flyem-cns-roi-7c971aa681da83f9a074a1f0e8ef60f4"
    return dict(
        title=title,
        dimensions={a: [8e-9, "m"] for a in "xyz"},
        position=[47470.59375, 54212.50390625, 66301.59375],
        crossSectionScale=1.1208,
        projectionOrientation=[0, 0.7071067690849304, 0.7071067690849304, 0],
        projectionScale=192860.74860894235,
        layers=[dict(type="image",
                     source=dict(url="precomputed://gs://cns-full-clahe",
                                 subsources=dict(default=True), enableDefaultSubsources=False),
                     tab="rendering", name="em"),
                shell(f"{roi}/brain-shell-smooth-linear", "brain-shell", "rendering"),
                shell(f"{roi}/vnc-shell", "vnc-shell", "segments")],
        showAxisLines=False, showSlices=False, prefetch=False, layout="3d")


def seg_layer(ids, name, colour):
    ids = [str(i) for i in ids]
    return dict(
        type="segmentation",
        source=[dict(url="dvid://https://emdata6-novran.janelia.org/4b2087c0fbe046bfaf0d60bc970e3e5d/"
                         "segmentation?dvid-service=https://ngsupport-bmcp5imp6q-uk.a.run.app",
                     subsources=dict(default=True, meshes=True), enableDefaultSubsources=False),
                "precomputed://https://ngsupport-bmcp5imp6q-uk.a.run.app/neuronjson_segment_properties/"
                "emdata6-novran.janelia.org/4b2087c0fbe046bfaf0d60bc970e3e5d/"
                "segmentation_annotations/type/group"],
        toolBindings=dict(Q="selectSegments"), tab="segments",
        segments=ids, segmentQuery=" ".join(ids),
        segmentColors={i: colour for i in ids}, name=name)


def point_layer(df, name, colour):
    return dict(
        type="annotation",
        source=dict(url="local://annotations",
                    transform=dict(outputDimensions={a: [8e-9, "m"] for a in "xyz"})),
        tool="annotatePoint",
        annotations=[dict(point=[float(r.x), float(r.y), float(r.z)], type="point", id=str(i))
                     for i, r in enumerate(df.itertuples())],
        annotationColor=colour,
        shader="void main() { setColor(defaultColor()); setPointMarkerSize(8.0); }",
        tab="annotations", name=f"{name} ({len(df)} syn)")


def axoaxo(pre, post, name, min_sites=20, open_scene=True):
    c = client()
    pre_ids, post_ids = resolve_ids(pre, c), resolve_ids(post, c)
    meta = neu.fetch_neurons(neu.NeuronCriteria(bodyId=pre_ids + post_ids), client=c)[0]
    tmap = dict(zip(meta.bodyId, meta.type))
    type_of = lambda b: tmap.get(b) if isinstance(tmap.get(b), str) and tmap.get(b) else "untyped"
    out = AXO_DIR / "scenes_navis" / name
    out.mkdir(parents=True, exist_ok=True)

    # every pre -> post synapse, located at the postsynaptic site
    conns = neu.fetch_synapse_connections(neu.NeuronCriteria(bodyId=pre_ids),
                                          neu.NeuronCriteria(bodyId=post_ids), client=c)
    rows = []
    for b in post_ids:
        n = fetch_skeleton(b, c)
        s = conns[conns.bodyId_post == b].copy()
        sp = split(n, min_sites)
        s["compartment"], s["seg_index"] = pd.NA, pd.NA
        if sp is not None and len(s):
            comp, seg = sp
            _, idx = cKDTree(n.nodes[["x", "y", "z"]].values).query(s[["x_post", "y_post", "z_post"]].values)
            s["compartment"] = comp.reindex(n.nodes.node_id.values[idx]).values
            s["seg_index"] = seg
        rows.append(s)
    syn = pd.concat(rows, ignore_index=True).rename(columns={
        "bodyId_pre": "pre", "bodyId_post": "post", "x_post": "x", "y_post": "y", "z_post": "z"})
    syn["domain"] = syn.compartment.map(DOMAIN_OF)
    syn["pre_type"] = syn.pre.map(type_of)
    syn["post_type"] = syn.post.map(type_of)
    syn = syn[["pre", "pre_type", "post", "post_type", "compartment", "domain", "seg_index", "x", "y", "z"]]

    summ = (syn.groupby(["post", "post_type", "pre_type"])
               .agg(n_syn=("domain", "size"),
                    n_axo_dendritic=("domain", lambda d: (d == "axo-dendritic").sum()),
                    n_axo_axonic=("domain", lambda d: (d == "axo-axonic").sum()),
                    n_axo_linker=("domain", lambda d: (d == "axo-linker").sum()),
                    seg_index=("seg_index", "first"))
               .reset_index())
    summ["pct_axo_axonic"] = 100 * summ.n_axo_axonic / summ.n_syn
    summ = summ.sort_values(["post_type", "post", "n_syn"], ascending=[True, True, False])
    unsplit = sorted(set(post_ids) - set(syn.post[syn.domain.notna()]))
    if unsplit:
        print(f"no split (< {min_sites} input or output sites, or no synapses from pre): "
              + ", ".join(map(str, unsplit)))
    syn.to_csv(out / "synapses.csv", index=False)
    summ.to_csv(out / "summary.csv", index=False)

    # scene: one segmentation layer per type, one point layer per pre type x domain
    scene = mcns_scene(name)
    shown_pre = [b for b in pre_ids if b in set(syn.pre)]
    for ids, role, pal in ((post_ids, "post", POST_PAL), (shown_pre, "pre", PRE_PAL)):
        by_type = pd.Series(ids).groupby(pd.Series(ids).map(type_of))
        for i, (tp, bb) in enumerate(by_type):
            scene["layers"].append(seg_layer(bb, f"{role} {tp} ({len(bb)})", pal[i % len(pal)]))
    lab = syn[syn.domain.notna()]
    for tp in sorted(lab.pre_type.unique()):
        for dm, col in DOMAIN_COLS.items():
            d = lab[(lab.pre_type == tp) & (lab.domain == dm)]
            if len(d):
                scene["layers"].append(point_layer(d, f"{tp} {dm}", col))
    scene["selectedLayer"] = dict(visible=True, layer=scene["layers"][3]["name"])
    if len(syn):
        scene["position"] = syn[["x", "y", "z"]].mean().tolist()
        span = (syn[["x", "y", "z"]].max() - syn[["x", "y", "z"]].min()).max()
        scene["projectionScale"] = max(2500, round(span * 1.8))

    js = json.dumps(scene, separators=(",", ":"))
    (out / "scene.json").write_text(js)
    url = "https://clio-ng.janelia.org/#!" + urllib.parse.quote(js, safe="")
    (out / "scene_url.txt").write_text(url)
    # ~60k-character URL; a .webloc opens it on double-click in Finder
    link = out / f"{name}.webloc"
    link.write_text('<?xml version="1.0" encoding="UTF-8"?>\n<plist version="1.0"><dict>'
                    f'<key>URL</key><string>{url.replace("&", "&amp;")}</string></dict></plist>')
    if open_scene:
        subprocess.run(["open", str(link)])

    print(f"{name}: {len(syn)} synapses, {len(post_ids)} post bodies, "
          f"{100 * (lab.domain == 'axo-axonic').mean():.1f}% axo-axonic "
          f"({100 * (lab.domain == 'axo-linker').mean():.1f}% on linker)")
    print("written to", out)
    print("scene link:", link)
    return syn, summ


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pre", nargs="+")
    ap.add_argument("--post", nargs="+")
    ap.add_argument("--name")
    ap.add_argument("--min-sites", type=int, default=20)
    ap.add_argument("--no-open", action="store_true")
    ap.add_argument("--set-token", nargs="?", const="", metavar="TOKEN",
                    help="save a new neuprint token to .env (prompts if TOKEN is omitted)")
    a = ap.parse_args()
    if a.set_token is not None:
        set_token(a.set_token)
    elif not (a.pre and a.post and a.name):
        ap.error("--pre, --post and --name are required (or use --set-token)")
    else:
        axoaxo(a.pre, a.post, a.name, a.min_sites, not a.no_open)
