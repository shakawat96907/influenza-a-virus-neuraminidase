import math
from pymol import cmd

cmd.load("results/structure/4B7R_tetramer.pdb", "na_tetramer")
cmd.remove("solvent")
cmd.hide("everything")
cmd.show("cartoon", "na_tetramer")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)

# Vibrant 4-subunit coloring
cmd.set_color("sub_a_blue", [0.11, 0.42, 0.78])
cmd.set_color("sub_b_teal", [0.09, 0.65, 0.56])
cmd.set_color("sub_c_coral", [0.86, 0.37, 0.34])
cmd.set_color("sub_d_gold", [0.89, 0.66, 0.18])

cmd.color("sub_a_blue", "chain A")
cmd.color("sub_b_teal", "chain B")
cmd.color("sub_c_coral", "chain C")
cmd.color("sub_d_gold", "chain D")
cmd.set("cartoon_transparency", 0.15, "not chain A")

# Top Season-Matched WHO Vaccine-Divergent Head Residues (>=15%) as sticks + spheres (CPK)
cmd.select("vac_div_sites", "chain A and resi 188+200+241+248+264+270+314+321+369+386+389+397+416+432+449+452")
cmd.show("sticks", "vac_div_sites")
cmd.set("stick_radius", 0.25, "vac_div_sites")
cmd.show("spheres", "vac_div_sites and name CA")
cmd.set("sphere_scale", 0.45, "vac_div_sites")
cmd.color("brightorange", "vac_div_sites and elem C")
cmd.util.cnc("vac_div_sites")

# Conserved catalytic pocket in deep crimson
cmd.select("active_site", "chain A and resi 118+151+152+225+275+277+293+368+402")
cmd.show("sticks", "active_site")
cmd.set("stick_radius", 0.22, "active_site")
cmd.color("ruby", "active_site and elem C")
cmd.util.cnc("active_site")

# Non-overlapping radial callout labels outside the tetramer silhouette
cx, cy, cz = cmd.centerofmass("na_tetramer")
landmarks = [
    (188, "T188 (35.1%)", 16.5, (-3.0, 5.5, 0.0)),
    (200, "S200 (27.7%)", 16.5, (2.5, 5.5, 0.0)),
    (264, "I264 (34.2%)", 15.5, (6.5, 1.5, 0.0)),
    (386, "K386 (23.9%)", 16.5, (7.5, 1.8, 0.0)),
    (416, "N416 (25.1%)", 16.0, (5.5, -6.5, 0.0)),
    (449, "D449 (36.0%)", 17.5, (-1.0, -8.5, 0.0))
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
cmd.zoom("na_tetramer", 12)
cmd.ray(2600, 2000)
cmd.png("results/structure/pymol_fig10_raw.png", dpi=300)
cmd.quit()
