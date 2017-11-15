//
// Created by wenhan on 9/4/17.
//

#ifndef DBSK2D_ISHOCK_COMPUTATION_DBSK2D_LOAD_OSF_PROCESS_H
#define DBSK2D_ISHOCK_COMPUTATION_DBSK2D_LOAD_OSF_PROCESS_H


#include "bpro1/bpro1_process.h"
#include "bpro1/bpro1_parameters.h"
#include "dbsk2d/dbsk2d_ishock_graph_sptr.h"
#include "dbsk2d/dbsk2d_ishock_belm.h"
#include "dbsk2d/dbsk2d_ishock_bline.h"
#include "dbsk2d/dbsk2d_ishock_bpoint.h"
#include "dbsk2d/pro/dbsk2d_shock_storage_sptr.h"

#include <vcl_cstring.h>
#include <vcl_iostream.h>
#include <vcl_fstream.h>

class dbsk2d_load_osf_process : public bpro1_process
{
private:
    dbsk2d_shock_storage_sptr _output_shock;

    dbsk2d_ishock_graph_sptr _ishock_graph;
    vcl_vector<dbsk2d_ishock_belm* > _belms;

    // get ishock node from node index
    dbsk2d_ishock_node* get_ishock_node_from_index(int index);

    // get bnd elm from bnd index
    dbsk2d_ishock_belm* get_bnd_from_index(int index);

    // get ishock edge from edge index
    dbsk2d_ishock_edge* get_ishock_edge_from_index(int index);


public:
    dbsk2d_load_osf_process();
    virtual ~dbsk2d_load_osf_process() {}
    virtual bpro1_process* clone() const;

    vcl_string name()
    {
        return "Load .osf file.";
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

    dbsk2d_shock_storage_sptr get_output_shock()
    {
        return _output_shock;
    }

    bool execute();

    bool finish()
    {
        vcl_cout << "Finish load .osf process!" << vcl_endl;
        return true;
    }

    bool load_osf(vcl_string file_name);

};


#endif //DBSK2D_ISHOCK_COMPUTATION_DBSK2D_LOAD_OSF_PROCESS_H
