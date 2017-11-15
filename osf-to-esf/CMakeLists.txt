cmake_minimum_required(VERSION 3.5)
project(osf_to_esf)


FIND_PACKAGE( VXL )

IF( VXL_FOUND )
    INCLUDE( ${VXL_CMAKE_DIR}/UseVXL.cmake )
ENDIF( VXL_FOUND )

IF( VXL_GEL_FOUND )
    INCLUDE_DIRECTORIES( ${VXL_GEL_INCLUDE_DIR} )
    MESSAGE("GEL FOUND")
ENDIF( VXL_GEL_FOUND )

IF( VXL_MUL_FOUND )
    INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )
    MESSAGE("MUL FOUND")
ENDIF( VXL_MUL_FOUND )

IF(VXL_BRL_FOUND )
    INCLUDE_DIRECTORIES( ${VXL_BRL_INCLUDE_DIR} )
    MESSAGE("BRL FOUND")
ENDIF( VXL_BRL_FOUND )


SET(SOURCE_FILES

        osf_to_esf_main.cpp
        dbsk2d_compute_ishock_process.cpp
        dbsk2d_compute_ishock_process.h
        dbsk2d_load_osf_process.cpp
        dbsk2d_load_osf_process.h

        bpro1/bpro1_parameters.cxx
        bpro1/bpro1_parameters.h
        bpro1/bpro1_process.cxx
        bpro1/bpro1_process.h
        bpro1/bpro1_storage_sptr.h
        bpro1/bpro1_parameters_sptr.h
        bpro1/bpro1_storage.cxx
        bpro1/bpro1_storage.h
        bpro1/bpro1_parameters.txx
        bpro1/bpro1_process_sptr.h

        vidpro1/storage/vidpro1_image_storage.cxx
        vidpro1/storage/vidpro1_image_storage.h
        vidpro1/storage/vidpro1_image_storage_sptr.h
        vidpro1/storage/vidpro1_vsol2D_storage.cxx
        vidpro1/storage/vidpro1_vsol2D_storage.h
        vidpro1/storage/vidpro1_vsol2D_storage_sptr.h

        vidpro1/pro/vidpro1_load_cem_process.cxx
        vidpro1/pro/vidpro1_load_cem_process.h

        dbgl/dbgl_eulerspiral_base.cxx
        dbgl/dbgl_eulerspiral_base.h
        dbgl/dbgl_arc.cxx
        dbgl/dbgl_arc.h
        dbgl/dbgl_eno_curve.cxx
        dbgl/dbgl_eno_curve.h
        dbgl/dbgl_poly_curve_line.cxx
        dbgl/dbgl_poly_curve_line.h
        dbgl/dbgl_arc_base.cxx
        dbgl/dbgl_arc_base.h
        dbgl/dbgl_poly_curve.cxx
        dbgl/dbgl_poly_curve_line.h

        dbgl/dbgl_fuzzy_boolean.h
        dbgl/dbgl_param_curve.cxx
        dbgl/dbgl_param_curve.h
        dbgl/algo/dbgl_closest_point.cxx
        dbgl/algo/dbgl_closest_point.h
        dbgl/algo/dbgl_biarc.cxx
        dbgl/algo/dbgl_biarc.h
        dbgl/algo/dbgl_circ_arc.cxx
        dbgl/algo/dbgl_circ_arc.h
        dbgl/algo/dbgl_distance.cxx
        dbgl/algo/dbgl_distance.h
        dbgl/algo/dbgl_intersect.cxx
        dbgl/algo/dbgl_intersect.h
        dbgl/algo/dbgl_eulerspiral.cxx
        dbgl/algo/dbgl_eulerspiral.h

        db3p/ray_tri/ray_tri.cxx
        db3p/ray_tri/ray_tri.h

        dbgrl/dbgrl_edge.txx
        dbgrl/dbgrl_edge.h
        dbgrl/dbgrl_graph.h
        dbgrl/dbgrl_graph.txx
        dbgrl/dbgrl_vertex.h
        dbgrl/dbgrl_vertex.txx

        dbnl/dbnl_solve_quadratic_equation.cxx
        dbnl/dbnl_solve_quadratic_equation.h
        dbnl/dbnl_solve_trigonometric_equation.cxx
        dbnl/dbnl_solve_trigonometric_equation.h
        dbnl/dbnl_fresnel.cxx
        dbnl/dbnl_fresnel.h

        dbnl/algo/dbnl_eno.cxx
        dbnl/algo/dbnl_eno.h

        dbsk2d/algo/dbsk2d_bnd_preprocess.cxx
        dbsk2d/algo/dbsk2d_bnd_preprocess.h
        dbsk2d/algo/dbsk2d_ishock_detector.cxx
        dbsk2d/algo/dbsk2d_ishock_detector.h
        dbsk2d/algo/dbsk2d_prune_ishock.cxx
        dbsk2d/algo/dbsk2d_prune_ishock.h
        dbsk2d/algo/dbsk2d_lagrangian_ishock_detector.cxx
        dbsk2d/algo/dbsk2d_lagrangian_ishock_detector.h
        dbsk2d/algo/dbsk2d_ishock_compute_intersection.cxx
        dbsk2d/algo/dbsk2d_ishock_compute_intersection.h
        dbsk2d/algo/dbsk2d_ishock_intersection_data.h
        dbsk2d/algo/dbsk2d_ishock_transform.cxx
        dbsk2d/algo/dbsk2d_ishock_transform.h
        dbsk2d/algo/dbsk2d_ishock_transform_sptr.h
        dbsk2d/algo/dbsk2d_ishock_gap_detector.cxx
        dbsk2d/algo/dbsk2d_ishock_gap_detector.h
        dbsk2d/algo/dbsk2d_ishock_gap_transform.cxx
        dbsk2d/algo/dbsk2d_ishock_gap_transform.h
        dbsk2d/algo/dbsk2d_ishock_grouping_transform.cxx
        dbsk2d/algo/dbsk2d_ishock_grouping_transform.h
        dbsk2d/algo/dbsk2d_sample_ishock.cxx
        dbsk2d/algo/dbsk2d_sample_ishock.h
        dbsk2d/algo/dbsk2d_bnd_preprocess_arcs.cxx
        dbsk2d/algo/dbsk2d_bnd_preprocess_common.cxx
        dbsk2d/algo/dbsk2d_bnd_preprocess_lines.cxx
        dbsk2d/algo/dbsk2d_bnd_preprocess_lines_vs_arcs.cxx
        dbsk2d/algo/dbsk2d_bnd_preprocess_points.cxx
        dbsk2d/algo/dbsk2d_lagrangian_ishock_detector_init.cxx
        dbsk2d/algo/dbsk2d_lagrangian_ishock_detector_prop.cxx
        dbsk2d/algo/dbsk2d_shock_transforms.cxx
        dbsk2d/algo/dbsk2d_shock_transforms.h
        dbsk2d/algo/dbsk2d_xshock_graph_fileio.cxx
        dbsk2d/algo/dbsk2d_xshock_graph_fileio.h

        dbsol/dbsol_interp_curve_2d.cxx
        dbsol/dbsol_interp_curve_2d.h
        dbsol/dbsol_interp_curve_2d_sptr.h
        dbsol/dbsol_file_io.cxx
        dbsol/dbsol_file_io.h
        dbsol/dbsol_curve_algs.cxx
        dbsol/dbsol_curve_algs.h

        dbsk2d/pro/dbsk2d_shock_storage.cxx
        dbsk2d/pro/dbsk2d_shock_storage.h
        dbsk2d/pro/dbsk2d_shock_storage_sptr.h
        dbsk2d/pro/dbsk2d_gap_transform_process.cxx
        dbsk2d/pro/dbsk2d_gap_transform_process.h
        dbsk2d/pro/dbsk2d_sample_ishock_process.cxx
        dbsk2d/pro/dbsk2d_sample_ishock_process.h
        dbsk2d/pro/dbsk2d_save_esf_process.cxx
        dbsk2d/pro/dbsk2d_save_esf_process.h

        dbsk2d/dbsk2d_boundary.cxx
        dbsk2d/dbsk2d_boundary.h
        dbsk2d/dbsk2d_boundary_sptr.h
        dbsk2d/dbsk2d_ishock_graph_sptr.h
        dbsk2d/dbsk2d_shock_graph.cxx
        dbsk2d/dbsk2d_shock_graph.h
        dbsk2d/dbsk2d_shock_graph_sptr.h
        dbsk2d/dbsk2d_bnd_utils.cxx
        dbsk2d/dbsk2d_bnd_utils.h
        dbsk2d/dbsk2d_ishock_arcarc.cxx
        dbsk2d/dbsk2d_ishock_arcarc.h
        dbsk2d/dbsk2d_ishock_arcarc_thirdorder.cxx
        dbsk2d/dbsk2d_ishock_arcarc_thirdorder.h
        dbsk2d/dbsk2d_ishock_barc.cxx
        dbsk2d/dbsk2d_ishock_barc.h
        dbsk2d/dbsk2d_ishock_bcurve.cxx
        dbsk2d/dbsk2d_ishock_bcurve.h
        dbsk2d/dbsk2d_ishock_belm.cxx
        dbsk2d/dbsk2d_ishock_belm.h
        dbsk2d/dbsk2d_ishock_bline.cxx
        dbsk2d/dbsk2d_ishock_bline.h
        dbsk2d/dbsk2d_ishock_bnd_key.cxx
        dbsk2d/dbsk2d_ishock_bnd_key.h
        dbsk2d/dbsk2d_ishock_bpoint.cxx
        dbsk2d/dbsk2d_ishock_bpoint.h
        dbsk2d/dbsk2d_ishock_contact.cxx
        dbsk2d/dbsk2d_ishock_contact.h
        dbsk2d/dbsk2d_ishock_edge.cxx
        dbsk2d/dbsk2d_ishock_edge.h
        dbsk2d/dbsk2d_ishock_elm.h
        dbsk2d/dbsk2d_ishock_graph.cxx
        dbsk2d/dbsk2d_ishock_graph.h
        dbsk2d/dbsk2d_ishock_graph_sptr.h
        dbsk2d/dbsk2d_ishock_linearc.cxx
        dbsk2d/dbsk2d_ishock_linearc.h
        dbsk2d/dbsk2d_ishock_lineline.cxx
        dbsk2d/dbsk2d_ishock_lineline.h
        dbsk2d/dbsk2d_ishock_lineline_thirdorder.cxx
        dbsk2d/dbsk2d_ishock_lineline_thirdorder.h
        dbsk2d/dbsk2d_ishock_node.cxx
        dbsk2d/dbsk2d_ishock_node.h
        dbsk2d/dbsk2d_ishock_pointarc.cxx
        dbsk2d/dbsk2d_ishock_pointarc.h
        dbsk2d/dbsk2d_ishock_pointarc_thirdorder.cxx
        dbsk2d/dbsk2d_ishock_pointarc_thirdorder.h
        dbsk2d/dbsk2d_ishock_pointline.cxx
        dbsk2d/dbsk2d_ishock_pointline.h
        dbsk2d/dbsk2d_ishock_pointpoint.cxx
        dbsk2d/dbsk2d_ishock_pointpoint.h
        dbsk2d/dbsk2d_ishock_utils.h
        dbsk2d/dbsk2d_lagrangian_cell_bnd.cxx
        dbsk2d/dbsk2d_lagrangian_cell_bnd.h
        dbsk2d/dbsk2d_lagrangian_cell_bnd_sptr.h
        dbsk2d/dbsk2d_rich_map_sptr.h
        dbsk2d/dbsk2d_shock_edge.cxx
        dbsk2d/dbsk2d_shock_edge.h
        dbsk2d/dbsk2d_shock_edge_sptr.h
        dbsk2d/dbsk2d_shock_grouped_ishock_edge.cxx
        dbsk2d/dbsk2d_shock_grouped_ishock_edge.h
        dbsk2d/dbsk2d_shock_grouped_ishock_edge_sptr.h
        dbsk2d/dbsk2d_shock_ishock_node.cxx
        dbsk2d/dbsk2d_shock_ishock_node.h
        dbsk2d/dbsk2d_shock_ishock_node_sptr.h
        dbsk2d/dbsk2d_shock_node.cxx
        dbsk2d/dbsk2d_shock_node.h
        dbsk2d/dbsk2d_shock_node_sptr.h
        dbsk2d/dbsk2d_utils.cxx
        dbsk2d/dbsk2d_utils.h
        dbsk2d/dbsk2d_bnd_edge.cxx
        dbsk2d/dbsk2d_bnd_edge.h
        dbsk2d/dbsk2d_bnd_edge_sptr.h
        dbsk2d/dbsk2d_bnd_cell.cxx
        dbsk2d/dbsk2d_bnd_cell.h
        dbsk2d/dbsk2d_bnd_cell_sptr.h
        dbsk2d/dbsk2d_bnd_contour.cxx
        dbsk2d/dbsk2d_bnd_contour.h
        dbsk2d/dbsk2d_bnd_contour_sptr.h
        dbsk2d/dbsk2d_bnd_vertex.cxx
        dbsk2d/dbsk2d_bnd_vertex.h
        dbsk2d/dbsk2d_bnd_vertex_sptr.h
        dbsk2d/dbsk2d_defines.h
        dbsk2d/dbsk2d_distance.cxx
        dbsk2d/dbsk2d_distance.h
        dbsk2d/dbsk2d_geometry_utils.cxx
        dbsk2d/dbsk2d_geometry_utils.h
        dbsk2d/dbsk2d_assert.cxx
        dbsk2d/dbsk2d_assert.h
        dbsk2d/dbsk2d_base_gui_geometry.h
        dbsk2d/dbsk2d_exception.cxx
        dbsk2d/dbsk2d_exception.h
        dbsk2d/dbsk2d_fuzzy_boolean.h
        dbsk2d/dbsk2d_shock_fragment.cxx
        dbsk2d/dbsk2d_shock_fragment.h
        dbsk2d/dbsk2d_shock_fragment_sptr.h
        dbsk2d/dbsk2d_rich_map.cxx
        dbsk2d/dbsk2d_rich_map.h
        dbsk2d/dbsk2d_containment_node.cxx
        dbsk2d/dbsk2d_containment_node.h
        dbsk2d/dbsk2d_containment_node_sptr.h
        dbsk2d/dbsk2d_file_io.cxx
        dbsk2d/dbsk2d_file_io.h
        dbsk2d/dbsk2d_transform_manager.cxx
        dbsk2d/dbsk2d_transform_manager.h
        dbsk2d/dbsk2d_xshock_edge.cxx
        dbsk2d/dbsk2d_xshock_edge.h
        dbsk2d/dbsk2d_xshock_sample.h
        dbsk2d/dbsk2d_xshock_sample_sptr.h
        dbsk2d/dbsk2d_xshock_node.cxx
        dbsk2d/dbsk2d_xshock_node.h

        params/vox_compute_ishock_params.cxx
        params/vox_compute_ishock_params.h
        params/vox_compute_ishock_params_sptr.h

        params/vbl_smart_ptr+vox_compute_ishock_params-.cxx

        dborl/algo/dborl_algo_params.cxx
        dborl/algo/dborl_algo_params.h
        dborl/dborl_algo_params_base.cxx
        dborl/dborl_algo_params_base.h
        dborl/dborl_algo_params_base_sptr.h
        dborl/dborl_parameter.cxx
        dborl/dborl_parameter.h
        dborl/dborl_evaluation.cxx
        dborl/dborl_evaluation.h
        dborl/dborl_image_description_utils.cxx
        dborl/dborl_image_description_utils.h
        dborl/dborl_category_info.cxx
        dborl/dborl_category_info.h
        dborl/dborl_category_info_sptr.h
        dborl/dborl_exp_stat_sptr.h
        dborl/dborl_image_description_sptr.h
        dborl/dborl_image_bbox_description.cxx
        dborl/dborl_image_bbox_description.h
        dborl/dborl_image_bbox_description_sptr.h
        dborl/dborl_image_data_description_base.cxx
        dborl/dborl_image_data_description_base.h
        dborl/dborl_image_data_description_base_sptr.h
        dborl/dborl_image_description.cxx
        dborl/dborl_image_description.h
        dborl/dborl_image_mask_description.cxx
        dborl/dborl_image_mask_description.h
        dborl/dborl_image_mask_description_sptr.h
        dborl/dborl_image_polygon_description.cxx
        dborl/dborl_image_polygon_description.h
        dborl/dborl_image_polygon_description_sptr.h
        dborl/dborl_description_base.h
        dborl/dborl_description_base_sptr.h
        dborl/algo/dborl_utilities.cxx
        dborl/algo/dborl_utilities.h
        )

# SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -pg )
# SET(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -pg)

AUX_SOURCE_DIRECTORY(dbsk2d/Templates SOURCE_FILES)
AUX_SOURCE_DIRECTORY(bpro1/Templates SOURCE_FILES)
AUX_SOURCE_DIRECTORY(vidpro1/Templates SOURCE_FILES)
AUX_SOURCE_DIRECTORY(dbsol/Templates SOURCE_FILES)
AUX_SOURCE_DIRECTORY(dborl/Templates SOURCE_FILES)

ADD_EXECUTABLE(osf_to_esf ${SOURCE_FILES})

LINK_DIRECTORIES(${VXL_DIR}/lib)

TARGET_LINK_LIBRARIES( osf_to_esf vul vgl vsol vil vcl vnl vnl_algo vbl vtol vil_io vnl_io bsol bsta vil3d vnl_io bbas_pro bdgl)



