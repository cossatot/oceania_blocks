using Revise

using Oiler

using CSV
using JSON
using DataFrames:eachrow
using DataFrames, DataFramesMeta
using Setfield

using PyPlot

# options
geol_slip_rate_weight = 2.
save_results = true

# load data

# ocn
ocn_block_file = "/home/itchy/research/geodesy/global_block_comps/oceania_blocks/block_data/oceania_blocks.geojson"
ocn_fault_file = "/home/itchy/research/geodesy/global_block_comps/oceania_blocks/block_data/oceania_faults.geojson"
ocn_geol_slip_rate_file = "/home/itchy/research/geodesy/global_block_comps/oceania_blocks/block_data/oceania_geol_slip_rates.geojson"

# chn
chn_block_file = "/home/itchy/research/geodesy/global_block_comps/china/block_data/chn_blocks.geojson"
chn_fault_file = "/home/itchy/research/geodesy/global_block_comps/china/block_data/chn_faults.geojson"
chn_slip_rate_file = "/home/itchy/research/geodesy/global_block_comps/china/block_data/geol_slip_rate_pts.geojson"

# jpn
izu_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/izu_slab2.geojson"

# phl

phl_block_file = "/home/itchy/research/geodesy/global_block_comps/phil_blocks/block_data/phl_blocks.geojson"
phl_fault_file = "/home/itchy/research/geodesy/global_block_comps/phil_blocks/block_data/phl_faults.geojson"

# glo
glo_block_file = "/home/itchy/research/geodesy/global_block_comps/global_scale_plates/global_scale_plates.geojson"
glo_fault_file = "/home/itchy/research/geodesy/global_block_comps/global_scale_plates/global_scale_faults.geojson"
glo_slip_rates_file = "/home/itchy/research/geodesy/global_block_comps/global_scale_plates/global_scale_slip_rates.geojson"

# tris
png_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/png_slab2.geojson"
phi_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/phi_slab2.geojson"
man_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/man_slab2.geojson"
cot_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/cot_slab2.geojson"
sul_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/sul_slab2.geojson"
sol_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/sol_slab2.geojson"
sum_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/sum_slab2.geojson"
van_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/van_slab2.geojson"
izu_tris_file = "/home/itchy/research/geodesy/global_block_comps/subduction/sub_tri_meshes/izu_slab2.geojson"

# bounds files
ocn_bounds_file = "../block_data/oceania_bounds.geojson"

# geodesy
gsrm_vels_file = "/home/itchy/research/geodesy/gsrm/gps/gps_su.geojson"
hsu_vels_file = "/home/itchy/research/geodesy/global_block_comps/phil_blocks/geod/hsu_luzon_campaign_vels.geojson"


@info "joining blocks"
glo_blocks = Oiler.IO.gis_vec_file_to_df(glo_block_file; fid_drop=["ant"])
chn_blocks = Oiler.IO.gis_vec_file_to_df(chn_block_file)
ocn_blocks = Oiler.IO.gis_vec_file_to_df(ocn_block_file); 
                                         #fid_drop=[ "ocn010"])
phl_blocks = Oiler.IO.gis_vec_file_to_df(phl_block_file)

block_df = vcat(ocn_blocks, 
                glo_blocks, 
                phl_blocks,
                chn_blocks; 
                cols=:union)

println("n blocks: ", size(block_df, 1))

@info "culling blocks"
println("n blocks before ", size(block_df, 1))
bound_df = Oiler.IO.gis_vec_file_to_df(ocn_bounds_file)
block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df; epsg=102016)
println("n blocks after ", size(block_df, 1))

@info "doing faults"
fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
    ocn_fault_file,
    chn_fault_file,
    phl_fault_file,
    glo_fault_file;
    block_df=block_df,
    lsd_default=10.,
    e_default=10.,
    check_blocks=true,
    #fid_drop=["ocnf054", "ocnf060"],
    subset_in_bounds=true)

println("n faults: ", length(faults))
println("n fault vels: ", length(fault_vels))

@info "doing non-fault block boundaries"
@time non_fault_bounds = Oiler.IO.get_non_fault_block_bounds(block_df, faults)
bound_vels = vcat(map(Oiler.Boundaries.boundary_to_vels, non_fault_bounds)...)
println("n non-fault-bound vels: ", length(bound_vels))

@info "doing geologic slip rates"
ocn_slip_rate_df = Oiler.IO.gis_vec_file_to_df(ocn_geol_slip_rate_file)
glo_slip_rate_df = Oiler.IO.gis_vec_file_to_df(glo_slip_rates_file)

geol_slip_rate_df = vcat(ocn_slip_rate_df, glo_slip_rate_df)

geol_slip_rate_df, geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
    geol_slip_rate_df,
    fault_df;
    weight=geol_slip_rate_weight
    )

println("n geol slip rates: ", length(geol_slip_rate_vels))

@info "doing GNSS"
gsrm_vel_df = Oiler.IO.gis_vec_file_to_df(gsrm_vels_file)
@time gsrm_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gsrm_vel_df, block_df;
    fix="su", epsg=102016,
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station
)

hsu_vel_df = Oiler.IO.gis_vec_file_to_df(hsu_vels_file)
@time hsu_vels = Oiler.IO.make_vels_from_gnss_and_blocks(hsu_vel_df, block_df;
    fix="ITRF08", epsg=102016,
    ve=:ve, vn=:vn, ee=:ee, en=:en, name=:site
)

gnss_vels = vcat(
             gsrm_vels,
             hsu_vels,
             )

println("n gnss vels: ", length(gnss_vels))

@info "putting it all together"
vels = vcat(fault_vels,
            gnss_vels,
            geol_slip_rate_vels
            )


@info "doing tris"
sol_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(sol_tris_file))
sul_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(sul_tris_file))
sum_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(sum_tris_file))
png_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(png_tris_file))
van_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(van_tris_file))
cot_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(cot_tris_file))
phi_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(phi_tris_file))
izu_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(izu_tris_file))

tris = vcat(
            sol_tris,
            van_tris,
            #png_tris,
            sum_tris,
            cot_tris,
            phi_tris,
            #sul_tris,
            #izu_tris,
            )


println("n total vels: ", length(vels))
vel_groups = Oiler.group_vels_by_fix_mov(vels)
tri_distance_weight = 5.


@info "Solving"
@time results = Oiler.solve_block_invs_from_vel_groups(vel_groups,
            faults=faults,
            elastic_floor=1e-5,
            tris=tris,
            tri_distance_weight=tri_distance_weight,
            regularize_tris=true,
            tri_priors=false,
            predict_vels=true,
            pred_se=true,
            check_closures=false,
            check_nans=true,
            sparse_lhs=true,
            constraint_method="kkt_sym",
            factorization="lu")


Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 directory="../web_viewer", ref_pole="su")


map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults, tris)
rates_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, 
                                           geol_slip_rate_vels, 
                                           fault_df, results)
show()

if save_results == true
    Oiler.IO.write_tri_results_to_gj(tris, results,
                                     "../results/ocn_tris.geojson",
                                     name="ocn tri results")
    Oiler.IO.write_fault_results_to_gj(results,
                                       "../results/ocn_faults.geojson",
                                       name="ocn_faults")
    Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups;
                                       name="../results/ocn_gnss_results.csv")
end

println("done!")
