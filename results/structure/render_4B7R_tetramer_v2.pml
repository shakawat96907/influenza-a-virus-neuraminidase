load results/structure/4B7R_tetramer.pdb, na_tetramer
remove solvent
hide everything
show cartoon, na_tetramer
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
color marine, chain A
color gray85, not chain A
set cartoon_transparency, 0.38, not chain A

# Key head-domain variants as sticks + spheres with element coloring (R3.8)
select head_variants, chain A and resi 200+222+264+270+339+382+389+416+449+452+453+469
show sticks, head_variants
set stick_radius, 0.22, head_variants
show spheres, head_variants and name CA
set sphere_scale, 0.40, head_variants
color tv_orange, head_variants and elem C
util.cnc("head_variants")

# Conserved catalytic pocket residues as sticks with element coloring
select active_site, chain A and resi 118+151+152+225+275+277+293+368+402
show sticks, active_site
set stick_radius, 0.20, active_site
color forest, active_site and elem C
util.cnc("active_site")

# Non-overlapping labels on key structural landmarks
set label_size, 16
set label_font_id, 7
set label_color, black
set label_outline_color, white
label chain A and resi 200 and name CA, "N200"
label chain A and resi 264 and name CA, "V264"
label chain A and resi 339 and name CA, "S339"
label chain A and resi 389 and name CA, "I389"
label chain A and resi 416 and name CA, "D416"
label chain A and resi 453 and name CA, "V453 (Interface)"
label chain A and resi 469 and name CA, "K469"

bg_color white
set ray_opaque_background, 1
set ray_shadows, 0
set antialias, 2
orient na_tetramer
zoom na_tetramer, 6
ray 2600, 1700
png results/structure/pymol_8C_raw.png, dpi=300
