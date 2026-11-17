load results/structure/4B7R_tetramer.pdb, na_tetramer
remove solvent
hide everything
show cartoon, na_tetramer
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
color marine, chain A
color gray85, not chain A
set cartoon_transparency, 0.35, not chain A

# Key head-domain variant residues as sticks + spheres with element coloring (R3.8)
select head_variants, chain A and resi 200+222+264+270+339+382+389+416+449+452+453+469
show sticks, head_variants
set stick_radius, 0.22, head_variants
show spheres, head_variants and name CA
set sphere_scale, 0.42, head_variants
color tv_orange, head_variants and elem C
util.cnc("head_variants")

# Conserved catalytic pocket residues as sticks with element coloring
select active_site, chain A and resi 118+151+152+225+275+277+293+368+402
show sticks, active_site
set stick_radius, 0.20, active_site
color forest, active_site and elem C
util.cnc("active_site")

set label_size, 18
set label_font_id, 7
set label_color, black
set label_outline_color, white
set label_position, [1.8, 1.8, 2.0]
label head_variants and name CA, "%s%s" % (resn, resi)

bg_color white
set ray_opaque_background, 1
set ray_shadows, 0
set antialias, 2
orient na_tetramer
zoom na_tetramer, 2
ray 2400, 2000
png results/figures/Figure8C_PyMOL_3D_Tetramer_preview.png, dpi=300
