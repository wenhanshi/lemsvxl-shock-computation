//
// Created by wenhan on 8/8/17.
//

#include <clocale>
#include "dbsk2d_save_isf_process.h"

dbsk2d_save_isf_process::dbsk2d_save_isf_process() : bpro1_process()
{
    _boundary = NULL;
    _ishock_graph = NULL;
    if( !parameters()->add( "Output .isf file <file_name>", "-isfoutput", bpro1_filepath("ishock.isf", "*isf")))
        vcl_cerr << "ERROR: Adding parameters in " __FILE__ << vcl_endl;

}

bpro1_process *dbsk2d_save_isf_process::clone() const
{
    return new dbsk2d_save_isf_process(*this);
}

vcl_vector< vcl_string > dbsk2d_save_isf_process::get_input_type()
{
    vcl_vector< vcl_string > to_return;
    to_return.clear();
    return to_return;
}

vcl_vector< vcl_string > dbsk2d_save_isf_process::get_output_type()
{
    vcl_vector< vcl_string > to_return;
    to_return.push_back("isf");
    return to_return;
}

bool dbsk2d_save_isf_process::execute()
{
    bpro1_filepath output_path;

    if(!_boundary)
    {
        vcl_cout << "ERROR: empty boundary in save .isf!" << vcl_endl;
        return false;
    }
    if(!_ishock_graph)
    {
        vcl_cout << "ERROR: empty ishock_graph in save .isf!" << vcl_endl;
        return false;
    }

    parameters()->get_value("-isfoutput", output_path);
    vcl_string output_file = output_path.path;

    return saveISF(output_file);
}

bool dbsk2d_save_isf_process::saveISF(vcl_string file_name)
{
    vcl_ofstream outfp(file_name.c_str(), vcl_ios::out);

    if (!outfp){
        vcl_cerr << " Error opening file  " << file_name.c_str() << vcl_endl;
        return false;
    }

    outfp << "# ==============ISHOCK FILE==============" << vcl_endl;

    outfp << "# ==============BOUNDARY=================" << vcl_endl;
    outfp << "# Boundary Elements:" << vcl_endl;
    outfp << "# [ID] [TYPE] [x1, y1] [x2, y2] ... [xn, yn]" << vcl_endl;
    outfp << "[BOUNDARY_BEGIN]" << vcl_endl;
    _boundary->update_belm_list();
    vcl_vector< dbsk2d_ishock_belm* > belms = _boundary->belm_list();
    for (vcl_vector< dbsk2d_ishock_belm* >::iterator it = belms.begin(); it != belms.end(); it++)
    {
        outfp << (*it)->id() << " ";
        outfp << (*it)->type() << " ";
        vcl_vector<vgl_point_2d<double> > pts = (*it)->ex_pts();
        for ( vcl_vector<vgl_point_2d<double> >::iterator pit = pts.begin(); pit != pts.end(); pit++)
        {
            outfp << (*pit).x() << " " << (*pit).y() << " ";
        }
        outfp << vcl_endl;
    }

    outfp << "# ==============ISHOCK NODE==============" << vcl_endl;
    outfp << "# IShock Nodes:" << vcl_endl;
    outfp << "# [ID] [x, y] [CHILD_ISHOCK_NODE_ID] [BND_1_ID] [BND_2_ID] ..." << vcl_endl;
    outfp << "[ISHOCK_NODE_BEGIN]" << vcl_endl;
    ishock_node_list nodes = _ishock_graph->all_nodes();
    for (ishock_node_list_iter inl = nodes.begin(); inl != nodes.end() ; inl++)
    {
        outfp << (*inl)->id() << " ";
        outfp << (*inl)->origin().x() << " " << (*inl)->origin().y() << " ";
        if ((*inl)->cShock() && (*inl)->cShock()->cSNode())
        {
            outfp << (*inl)->cShock()->cSNode()->id() << " ";
        }
        else
        {
            outfp << "-1 ";
        }
        ishock_node_belm_list bnds = (*inl)->bndList();
        for (ishock_node_belm_list_iter inbl = bnds.begin(); inbl != bnds.end(); inbl++)
        {
            outfp << (*inbl).second->id() << " ";
        }
        outfp << vcl_endl;
    }

    outfp << "# ==============ISHOCK EDGE==============" << vcl_endl;
    outfp << "# IShock Edges:" << vcl_endl;
    outfp << "# [ID] [TYPE] [FROM_NODE_ID] [TO_NODE_ID] [LEFT_BND_ID] [RIGHT_BND_ID]" << vcl_endl;
    outfp << "[ISHOCK_EDGE_BEGIN]" << vcl_endl;
    ishock_edge_list edges = _ishock_graph->all_edges();
    for (ishock_edge_list_iter iel = edges.begin(); iel != edges.end(); iel++)
    {
        outfp << (*iel)->id() << " ";
        outfp << (*iel)->type() << " ";
        if ((*iel)->pSNode())
            outfp << (*iel)->pSNode()->id() << " ";
        else
            outfp << "-1 ";
        if ((*iel)->cSNode())
            outfp << (*iel)->cSNode()->id() << " ";
        else
            outfp << "-1 ";
        outfp << (*iel)->lBElement()->id() << " " << (*iel)->rBElement()->id() ;
        outfp << vcl_endl;
    }

    outfp.close();
    return true;
}

