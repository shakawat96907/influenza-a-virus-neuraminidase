import math
from pymol import cmd

cmd.load("results/structure/4B7R_tetramer.pdb", "na_tetramer")
cmd.remove("solvent")
cmd.hide("everything")
cmd.show("cartoon", "na_tetramer")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)

# Vibrant 4-subunit coloring for a lively biological homotetramer
cmd.color("marine", "chain A")
cmd.color("deepteal", "chain B")
cmd.color("mediumpurple", "chain C")
cmd.color("palecyan", "chain D")
cmd.set("cartoon_transparency", 0.18, "not chain A")

# Key head-domain variants as vibrant orange sticks + spheres with element coloring (R3.8)
cmd.select("head_variants", "chain A and resi 200+222+264+270+339+382+389+416+449+452+453+469")
cmd.show("sticks", "head_variants")
cmd.set("stick_radius", 0.25, "head_variants")
cmd.show("spheres", "head_variants and name CA")
cmd.set("sphere_scale", 0.45, "head_variants")
cmd.color("brightorange", "head_variants and elem C")
cmd.util.cnc("head_variants")

# Conserved catalytic active-site pocket in vibrant ruby red
cmd.select("active_site", "chain A and resi 118+151+152+225+275+277+293+368+402")
cmd.show("sticks", "active_site")
cmd.set("stick_radius", 0.22, "active_site")
cmd.color("ruby", "active_site and elem C")
cmd.util.cnc("active_site")

# Radial callout labels outside the structure with subtle tinted background
cx, cy, cz = cmd.centerofmass("chain A")
landmarks = [
    (200, "N200", 11.5, (0.0, 4.0, 0.0)),
    (264, "V264", 11.5, (4.0, -2.5, 0.0)),
    (339, "S339", 11.0, (2.5, 0.0, 0.0)),
    (416, "D416", 11.0, (3.0, 0.0, 0.0)),
    (453, "V453 (Interface)", 12.0, (0.0, -4.0, 0.0)),
    (469, "K469", 11.5, (-3.5, -2.5, 0.0))
]

cmd.set("label_size", 15)
cmd.set("label_font_id", 7)
cmd.set("label_color", "black")
cmd.set("label_bg_color", "white")
cmd.set("label_bg_transparency", 0.15)
cmd.set("label_bg_outline", 1)
cmd.set("dash_gap", 0.0)
cmd.set("dash_width", 1.6)
cmd.set("dash_color", "gray30")

for resi, text, dist_out, nudge in landmarks:
    coords = cmd.get_coords(f"chain A and resi {resi} and name CA", 1)
    if coords is not None and len(coords) > 0:
        x, y, z = float(coords), float(coords), float(coords)
        dx, dy, dz = x - cx, y - cy, z - cz
        norm = math.sqrt(dx*dx + dy*dy + dz*dz) or 1.0
        px = x + (dx / norm) * dist_out + nudge
        py = y + (dy / norm) * dist_out + nudge
        pz = z + (dz / norm) * dist_out + nudge
        pname = f"lbl_{resi}"
        dname = f"line_{resi}"
        cmd.pseudoatom(pname, pos=[px, py, pz])
        cmd.label(pname, f'"{text}"')
        cmd.distance(dname, f"chain A and resi {resi} and name CA", pname)
        cmd.hide("labels", dname)

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.orient("na_tetramer")
cmd.zoom("na_tetramer", 8)
cmd.ray(2600, 1700)
cmd.png("results/structure/pymol_8C_raw.png", dpi=300)
cmd.quit()
