//
// Created by wenhan on 8/14/17.
//

#ifndef DBSK2D_ISHOCK_COMPUTE_DBSK2D_LOAD_ISF_PROCESS_H
#define DBSK2D_ISHOCK_COMPUTE_DBSK2D_LOAD_ISF_PROCESS_H


#include "bpro1/bpro1_process.h"
#include "bpro1/bpro1_parameters.h"
#include "dbsk2d/dbsk2d_ishock_graph_sptr.h"
#include "dbsk2d/dbsk2d_ishock_belm.h"
#include "dbsk2d/dbsk2d_ishock_bline.h"
#include "dbsk2d/dbsk2d_ishock_bpoint.h"

#include <vcl_cstring.h>
#include <vcl_iostream.h>
#include <vcl_fstream.h>

class dbsk2d_load_isf_process : public bpro1_process
{
private:
    dbsk2d_ishock_graph_sptr _ishock_graph;
    vcl_vector<dbsk2d_ishock_belm* > _belms;

    // from line read, get its boundary info from this line
    void get_bnds(char buf[], vcl_vector<int>& bnd_ids);

    // get ishock node from node index
    dbsk2d_ishock_node* get_ishock_node_from_index(int index);

    // get ishock edge from edge index
    dbsk2d_ishock_belm* get_bnd_from_index(int index);

public:
    dbsk2d_load_isf_process();
    virtual ~dbsk2d_load_isf_process() {}
    virtual bpro1_process* clone() const;

    vcl_string name()
    {
        return "Load .isf file.";
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

    void set_ishock_graph(dbsk2d_ishock_graph_sptr ishock_graph)
    {
        _ishock_graph = ishock_graph;
    }

    dbsk2d_ishock_graph_sptr get_ishock_graph(void)
    {
        return _ishock_graph;
    }

    void set_belms(vcl_vector<dbsk2d_ishock_belm* >& belms)
    {
        _belms = belms;
    }

    vcl_vector<dbsk2d_ishock_belm* > get_belms()
    {
        return _belms;
    }

    bool execute();

    bool finish()
    {
        vcl_cout << "Finish load .isf process!" << vcl_endl;
        return true;
    }

    bool loadISF(vcl_string file_name);

};


#endif //DBSK2D_ISHOCK_COMPUTE_DBSK2D_LOAD_ISF_PROCESS_H
