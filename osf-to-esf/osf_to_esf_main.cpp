
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
#include "dbsk2d/pro/dbsk2d_gap_transform_process.h"
#include "dbsk2d/pro/dbsk2d_sample_ishock_process.h"
#include "dbsk2d/pro/dbsk2d_save_esf_process.h"
#include "dbsk2d_load_osf_process.h"

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

    // test .osf loader
    //******************** IShock Input *******************************
    // Add by Wenhan: 15/08/17
    // for testing ishock input
    vcl_cout << "************ IShock Input *************" << vcl_endl;
    dbsk2d_load_osf_process load_osf_pro;
    set_process_parameters_of_bpro1(*params,
                                    load_osf_pro,
                                    params->tag_load_osf_);
    if(load_osf_pro.execute())
        load_osf_pro.finish();
    else
        vcl_cout << "ERROR: LOAD OSF CRUSHED." << vcl_endl;

    dbsk2d_shock_storage_sptr output_shock = dbsk2d_shock_storage_new();
    output_shock = load_osf_pro.get_output_shock();
    output_shock->set_image(img_sptr);

    //********************   Gap Transform    *******************************

    // Perform Gap Transform on Intrinsinc Shock Graph
    vcl_vector<bpro1_storage_sptr> shock_gt_results;

    if ( params->gap_transform_())
    {

        vcl_cout<<"************  Gap Transform  *************"<<vcl_endl;

        dbsk2d_gap_transform_process gt_pro;
        set_process_parameters_of_bpro1(*params,
                                        gt_pro,
                                        params->tag_gap_transform_);


        // Before we start the process lets clean input output
        gt_pro.clear_input();
        gt_pro.clear_output();


        // Use input from ishock computation
        gt_pro.add_input(output_shock);
        bool gt_status = gt_pro.execute();
        gt_pro.finish();

        // Grab output from symbolic edge linking
        if ( gt_status )
        {
            shock_gt_results = gt_pro.get_output();
        }

        //Clean up after ourselves
        gt_pro.clear_input();
        gt_pro.clear_output();

        // It has two outputs
        if (shock_gt_results.size() != 2)
        {
            vcl_cerr << "Shock after gap transform is Invalid!"
                     << vcl_endl;
            return 1;
        }

    }

    //******************** Sample Shocks  ********************************
    vcl_cout<<"************  Sampling Shock *************"<<vcl_endl;

    dbsk2d_sample_ishock_process sample_sg_pro;
    set_process_parameters_of_bpro1(*params,
                                    sample_sg_pro,
                                    params->tag_sample_shock_);


    // Before we start the process lets clean input output
    sample_sg_pro.clear_input();
    sample_sg_pro.clear_output();

    // Use input from either ishock computation or gap_transform
    if ( params->gap_transform_())
    {
        sample_sg_pro.add_input(shock_gt_results[0]);
    }
    else
    {
        sample_sg_pro.add_input(output_shock);
    }

    // Kick of process
    bool sample_status = sample_sg_pro.execute();
    sample_sg_pro.finish();

    // Grab output from sampling
    vcl_vector<bpro1_storage_sptr> sample_shock_results;
    if ( sample_status )
    {
        sample_shock_results   = sample_sg_pro.get_output();
    }

    //Clean up after ourselves
    sample_sg_pro.clear_input();
    sample_sg_pro.clear_output();

    if (sample_shock_results.size() != 1)
    {
        vcl_cerr << "Sampling of Intrinsinc Shock Computation Failed"
                 << vcl_endl;
        return 1;
    }


    //******************** Save Shocks   ********************************
    vcl_cout<<"************   Saving  Shock *************"<<vcl_endl;

    dbsk2d_save_esf_process save_sg_pro;

    vcl_string output_file;
    if (params->save_to_object_folder_())
    {
        output_file = params->output_shock_folder_() + "/";
    }
    else
    {
        output_file = params->input_object_dir_() + "/";
    }

    if (!vul_file::exists(output_file))
    {
        vul_file::make_directory(output_file);

    }

    output_file = output_file + params->input_object_name_()+
                  params->output_extension_();

    bpro1_filepath output(output_file,params->output_extension_());

    save_sg_pro.parameters()->set_value("-esfoutput",output);

    // Before we start lets clean input output
    save_sg_pro.clear_input();
    save_sg_pro.clear_output();

    save_sg_pro.add_input(sample_shock_results[0]);
    bool status = save_sg_pro.execute();
    save_sg_pro.finish();

    //Clean up after ourselves
    save_sg_pro.clear_input();
    save_sg_pro.clear_output();

    if ( !status )
    {
        vcl_cerr << "Problems in saving extrinsinc shock file: "
                 << output_file << vcl_endl;
        return 1;

    }

    double vox_time = t.real()/1000.0;
    t.mark();
    vcl_cout<<vcl_endl;
    vcl_cout<<"************ Time taken: "<<vox_time<<" sec"<<vcl_endl;

    return 0;


}