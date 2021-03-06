
// By Wenhan 26/07

#include "params/vox_compute_ishock_params.h"
#include "params/vox_compute_ishock_params_sptr.h"
#include "dborl/algo/dborl_utilities.h"
#include <vul/vul_file.h>
#include <vul/vul_timer.h>
#include <vul/vul_file_iterator.h>
#include <bbas/bsol/bsol_algs.h>
#include <vsol/vsol_polygon_2d.h>
#include <vsol/vsol_box_2d.h>
#include <vil/vil_load.h>

#include "vidpro1/pro/vidpro1_load_cem_process.h"
#include "vidpro1/storage/vidpro1_vsol2D_storage_sptr.h"
#include "vidpro1/storage/vidpro1_vsol2D_storage.h"
#include "vidpro1/storage/vidpro1_image_storage_sptr.h"
#include "vidpro1/storage/vidpro1_image_storage.h"

#include "dbsk2d/pro/dbsk2d_shock_storage.h"
#include "dbsk2d/pro/dbsk2d_shock_storage_sptr.h"

#include "dbsk2d_compute_ishock_process.h"
#include "dbsk2d_save_isf_process.h"
#include "dbsk2d_load_isf_process.h"
#include "dbsk2d_save_osf_process.h"
//#include "dbsk2d/pro/dbsk2d_gap_transform_process.h"
//#include "dbsk2d/pro/dbsk2d_sample_ishock_process.h"
//#include "dbsk2d/pro/dbsk2d_save_esf_process.h"

int main(int argc, char *argv[]) {
    // Let time how long this takes
    // Start timer
    vul_timer t;

    // construct parameters with the default values;
    vox_compute_ishock_params_sptr params =
            new vox_compute_ishock_params("dborl_compute_ishock");

    // parse the command line arguments
    if (!params->parse_command_line_args(argc, argv))
        return 1;

    //: always print the params file if an executable to work with ORL web
    // interface
    if (!params->print_params_xml(params->print_params_file())) {
        vcl_cerr << "problems in writing params file to: "
                 << params->print_params_file() << vcl_endl;
    }

    // exit if there is nothing else to do
    if (params->exit_with_no_processing() || params->print_params_only()) {
        return 0;
    }

    //: always call this method to actually parse the input parameter file
    // whose name is extracted from the command line
    if (!params->parse_input_xml()) {
        return 1;
    }

    //Determine which input object we are going to use
    //Either from the input_object_dir or the associated file
    //The associated file always takes precendence
    vcl_string input_vsol_fn;


        // Use the database entries
    input_vsol_fn = params->input_object_dir_() + "/" + params->input_object_name_() + params->input_contour_extension_();


    if (!vul_file::exists(input_vsol_fn)) {
        vcl_cerr << "Cannot find contour file: " << input_vsol_fn << vcl_endl;
        return 1;
    }

    vcl_string input_contour_extension = vul_file::extension(input_vsol_fn);

    // Create output storage
    vcl_vector <bpro1_storage_sptr> vsol_contour;

    //Call appropriate process to load file
    if (input_contour_extension == ".cem" ||
        input_contour_extension == ".cemv") {

        // Call vidpro load cem process
        vidpro1_load_cem_process load_pro_cem;

        bpro1_filepath input(input_vsol_fn, input_contour_extension);

        load_pro_cem.parameters()->set_value("-ceminput", input);

        // Before we start the process lets clean input output
        load_pro_cem.clear_input();
        load_pro_cem.clear_output();

        // Pass in input vsol string
        bool load_cem_status = load_pro_cem.execute();
        load_pro_cem.finish();

        // Grab output from symbolic edge linking
        if (load_cem_status) {
            vsol_contour = load_pro_cem.get_output();
        }

        //Clean up after ourselves
        load_pro_cem.clear_input();
        load_pro_cem.clear_output();

        vcl_cout << input_vsol_fn << " loaded." << vcl_endl;
    } else {
        vcl_cerr << "Unknown input type: " <<
                 input_contour_extension << " Quit now" << vcl_endl;
        return 1;
    }

    if (vsol_contour.size() != 1) {
        vcl_cerr << "Problem loading file, " << input_vsol_fn <<
                 "Could not get vsol data structure" << vcl_endl;
        return 1;
    }

    // Grab the underlying contours
    vidpro1_vsol2D_storage_sptr vsol_contour_storage =
            vidpro1_vsol2D_storage_new();
    vsol_contour_storage.vertical_cast(vsol_contour[0]);

    if (params->add_bbox_()) {
        vcl_cout << "************  Compute  Bbox  *************" << vcl_endl;

        // Grab the underlying contours
        vidpro1_vsol2D_storage_sptr vsol_contour_storage;
        vsol_contour_storage.vertical_cast(vsol_contour[0]);

        // create new bounding box
        vsol_box_2d_sptr bbox = new vsol_box_2d();

        // Determine largest bounding box
        vcl_vector <vsol_spatial_object_2d_sptr> vsol_list =
                vsol_contour_storage->all_data();

        for (unsigned int b = 0; b < vsol_list.size(); b++) {
            bbox->grow_minmax_bounds(vsol_list[b]->get_bounding_box());
        }

        // Enlarge bounding box from size
        // Calculate xcenter, ycenter
        double xcenter = bbox->width() / 2.0;
        double ycenter = bbox->height() / 2.0;

        // Translate to center and scale
        double xmin_scaled = ((bbox->get_min_x() - xcenter) * 2) + xcenter;
        double ymin_scaled = ((bbox->get_min_y() - ycenter) * 2) + ycenter;
        double xmax_scaled = ((bbox->get_max_x() - xcenter) * 2) + xcenter;
        double ymax_scaled = ((bbox->get_max_y() - ycenter) * 2) + ycenter;

        bbox->add_point(xmin_scaled, ymin_scaled);
        bbox->add_point(xmax_scaled, ymax_scaled);

        vcl_cout << "bbox (minx, miny) (maxx, maxy) (width, height): "
                 << "(" << bbox->get_min_x() << ", " << bbox->get_min_y()
                 << ") (" << bbox->get_max_x() << ", " << bbox->get_max_y()
                 << ") ("
                 << bbox->width() << ", "
                 << bbox->height() << ")" << vcl_endl;

        // Add to vidpro storage this new bounding box
        vsol_polygon_2d_sptr box_poly = bsol_algs::poly_from_box(bbox);
        vsol_contour_storage->add_object(box_poly->cast_to_spatial_object());
    }

    //Make sure the input image exists
    vcl_string input_image_fn = params->input_object_dir_() + "/"
                                + params->input_object_name_() + params->input_image_extension_();

    if (!vul_file::exists(input_image_fn))
    {
        vcl_cerr << "Cannot find image file: " << input_image_fn << vcl_endl;
        return 1;
    }

    // Grab image
    vil_image_resource_sptr img_sptr =
            vil_load_image_resource(input_image_fn.c_str());
    if (!img_sptr)
    {
        vcl_cerr << "Cannot load image: " << input_image_fn << vcl_endl;
        return 1;
    }

    // Create vid pro storage
    vidpro1_image_storage_sptr inp = new vidpro1_image_storage();
    inp->set_image(img_sptr);

    vcl_cout << input_image_fn << " loaded." << vcl_endl;


    //******************** IShock Computation *******************************
    vcl_cout << "************ Computing Shock *************" << vcl_endl;

    dbsk2d_compute_ishock_process shock_pro;
    set_process_parameters_of_bpro1(*params,
                                    shock_pro,
                                    params->tag_compute_ishock_);

    // Before we start the process lets clean input output
    shock_pro.clear_input();
    shock_pro.clear_output();

    // Use input from edge detection
    shock_pro.add_input(inp);
    shock_pro.add_input(vsol_contour_storage);
    bool ishock_status = shock_pro.execute();
    shock_pro.finish();

    // Grab output from symbolic edge linking
    vcl_vector <bpro1_storage_sptr> shock_results;

    // If ishock status is bad we will keep iterating with noise till we get
    // a valid shock computation otherwise call it quits
    if (!ishock_status) {
        // Add noise to parameter set
        shock_pro.parameters()->set_value("-b_noise", true);

        // Clean up before we start running
        shock_pro.clear_input();
        shock_pro.clear_output();

        unsigned int i(0);
        unsigned int num_iterations = params->num_iter_();

        for (; i < num_iterations; ++i) {
            vcl_cout << vcl_endl;
            vcl_cout << "************ Retry Compute Shock,iter: "
                     << i + 1 << " *************" << vcl_endl;

            // Add inputs : vosl and image
            shock_pro.add_input(inp);
            shock_pro.add_input(vsol_contour_storage);

            // Kick off process again
            ishock_status = shock_pro.execute();
            shock_pro.finish();

            if (ishock_status) {
                // We have produced valid shocks lets quit
                break;
            }

            // Clean up after ourselves
            shock_pro.clear_input();
            shock_pro.clear_output();
        }
    }

    if (ishock_status) {
        shock_results = shock_pro.get_output();

        // Clean up after ourselves
        shock_pro.clear_input();
        shock_pro.clear_output();
    }


    //******************** IShock Output *******************************
    // Add by Wenhan: 09/08/17
    // for ishock output (.isf)
    vcl_cout << "************ IShock Output *************" << vcl_endl;
    dbsk2d_save_isf_process save_isf_pro;
    set_process_parameters_of_bpro1(*params,
                                    save_isf_pro,
                                    params->tag_save_isf_);

    dbsk2d_shock_storage_sptr ishock_sptr;
    ishock_sptr.vertical_cast(shock_results[0]);

    save_isf_pro.set_boundary(ishock_sptr->get_boundary());
    save_isf_pro.set_ishock_graph(ishock_sptr->get_ishock_graph());

    if(save_isf_pro.execute())
        save_isf_pro.finish();
    else
        vcl_cout << "ERROR: SAVE ISF CRUSHED." << vcl_endl;

    // test .isf loader
    //******************** IShock Input *******************************
    // Add by Wenhan: 15/08/17
    // for testing ishock input
    vcl_cout << "************ IShock Input *************" << vcl_endl;
    dbsk2d_load_isf_process load_isf_pro;
    set_process_parameters_of_bpro1(*params,
                                    load_isf_pro,
                                    params->tag_load_isf_);
    if(load_isf_pro.execute())
        load_isf_pro.finish();
    else
        vcl_cout << "ERROR: LOAD ISF CRUSHED." << vcl_endl;

    //******************** IShock Output (OSF) *******************************
    // Add by Wenhan: 04/09/17
    // for ishock output (.osf)
    vcl_cout << "************ IShock Output *************" << vcl_endl;
    dbsk2d_save_osf_process save_osf_pro;
    set_process_parameters_of_bpro1(*params,
                                    save_osf_pro,
                                    params->tag_save_osf_);

    save_osf_pro.set_boundary(ishock_sptr->get_boundary());
    save_osf_pro.set_ishock_graph(ishock_sptr->get_ishock_graph());

    if(save_osf_pro.execute())
        save_osf_pro.finish();
    else
        vcl_cout << "ERROR: SAVE OSF CRUSHED." << vcl_endl;

}