import math
from pymol import cmd

cmd.load("results/structure/4B7R_tetramer.pdb", "na_tetramer")
cmd.remove("solvent")
cmd.hide("everything")
cmd.show("cartoon", "na_tetramer")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)

# Rich, vivid 4-subunit journal palette (no dull gray/pale tones)
cmd.set_color("sub_a_blue", [0.11, 0.42, 0.78])
cmd.set_color("sub_b_teal", [0.09, 0.65, 0.56])
cmd.set_color("sub_c_coral", [0.86, 0.37, 0.34])
cmd.set_color("sub_d_gold", [0.89, 0.66, 0.18])

cmd.color("sub_a_blue", "chain A")
cmd.color("sub_b_teal", "chain B")
cmd.color("sub_c_coral", "chain C")
cmd.color("sub_d_gold", "chain D")
cmd.set("cartoon_transparency", 0.15, "not chain A")

# Key head-domain variants as bright orange-red sticks + spheres with element coloring (R3.8)
cmd.select("head_variants", "chain A and resi 200+222+264+270+339+382+389+416+449+452+453+469")
cmd.show("sticks", "head_variants")
cmd.set("stick_radius", 0.26, "head_variants")
cmd.show("spheres", "head_variants and name CA")
cmd.set("sphere_scale", 0.46, "head_variants")
cmd.color("brightorange", "head_variants and elem C")
cmd.util.cnc("head_variants")

# Conserved catalytic pocket in deep crimson
cmd.select("active_site", "chain A and resi 118+151+152+225+275+277+293+368+402")
cmd.show("sticks", "active_site")
cmd.set("stick_radius", 0.24, "active_site")
cmd.color("crimson", "active_site and elem C")
cmd.util.cnc("active_site")

# Push all callout labels cleanly outside the tetramer silhouette
cx, cy, cz = cmd.centerofmass("na_tetramer")
landmarks = [
    (200, "N200", 16.0, (2.0, 5.0, 0.0)),
    (264, "V264", 15.0, (6.0, 1.5, 0.0)),
    (339, "S339", 15.0, (6.5, -1.0, 0.0)),
    (416, "D416", 15.0, (5.5, -4.5, 0.0)),
    (453, "V453 (Interface)", 17.5, (4.0, -8.5, 0.0)),
    (469, "K469", 16.0, (6.5, 4.5, 0.0))
]

cmd.set("label_size", 15)
cmd.set("label_font_id", 7)
cmd.set("label_color", "black")
cmd.set("label_bg_color", "white")
cmd.set("label_bg_transparency", 0.05)
cmd.set("label_bg_outline", 1)
cmd.set("dash_gap", 0.0)
cmd.set("dash_width", 1.8)
cmd.set("dash_color", "gray25")

for resi, text, dist_out, nudge in landmarks:
    coords = cmd.get_coords(f"chain A and resi {resi} and name CA", 1)
    if coords is not None and len(coords) > 0:
        x, y, z = coords.flatten().tolist()
        dx, dy, dz = x - cx, y - cy, z - cz
        norm = math.sqrt(dx*dx + dy*dy + dz*dz) or 1.0
        nx, ny, nz = nudge
        px = x + (dx / norm) * dist_out + nx
        py = y + (dy / norm) * dist_out + ny
        pz = z + (dz / norm) * dist_out + nz
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
cmd.zoom("na_tetramer", 11)
cmd.ray(2600, 1800)
cmd.png("results/structure/pymol_8C_raw.png", dpi=300)
cmd.quit()
