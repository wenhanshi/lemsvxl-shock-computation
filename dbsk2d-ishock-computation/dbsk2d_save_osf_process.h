//
// Created by wenhan on 9/4/17.
//

#ifndef DBSK2D_ISHOCK_COMPUTATION_DBSK2D_SAVE_OSF_PROCESS_H
#define DBSK2D_ISHOCK_COMPUTATION_DBSK2D_SAVE_OSF_PROCESS_H

#include "bpro1/bpro1_process.h"
#include "bpro1/bpro1_parameters.h"
#include "dbsk2d/dbsk2d_ishock_belm.h"
#include "dbsk2d/dbsk2d_ishock_bline.h"
#include "dbsk2d/dbsk2d_boundary.h"
#include <vcl_iomanip.h>
#include <vcl_vector.h>
#include <vsl/vsl_binary_io.h>

class dbsk2d_save_osf_process : public bpro1_process
{
private:
    // we need recreate output_shock structure
    dbsk2d_boundary_sptr _boundary;
    dbsk2d_ishock_graph_sptr _ishock_graph;

public:

    dbsk2d_save_osf_process();

    virtual ~dbsk2d_save_osf_process() {}

    vcl_string name()
    {
        return "Save .osf file.";
    }

    virtual bpro1_process* clone() const;

    void set_boundary(dbsk2d_boundary_sptr boundary)
    {
        _boundary = boundary;
    }

    void set_ishock_graph(dbsk2d_ishock_graph_sptr ishock_graph)
    {
        _ishock_graph = ishock_graph;
    }

    vcl_vector< vcl_string > get_input_type();
    vcl_vector< vcl_string > get_output_type();

    int input_frames()
    {
        return 1;
    }
    int output_frames()
    {
        return 1;
    }

    bool execute();

    bool finish()
    {
        vcl_cout << "Finish save .osf process!" << vcl_endl;
    }

    bool save_osf(vcl_string file_name);
};


#endif //DBSK2D_ISHOCK_COMPUTATION_DBSK2D_SAVE_OSF_PROCESS_H
