//
// Created by wenhan on 8/14/17.
//

/*
 * For this loader, we modified some func and private members in key classes, to make the .isf easier to use.
 * We provide three groups of infomation in .isf file:
 * 1. boundaries elements
 * 2. ishock nodes
 * 3. ishock edges
 *
 * In boundaries, due to line fitting, we only read following types of boundaries:
 * - bpoint
 * - bline
 * with the message:
 * [ID] [TYPE] [x1, y1] [x2, y2] ... [xn, yn]
 * TODO: add barc in .isf saver and loader, maybe use arc fitting
 *
 * In ishock node, the shock type must be 1 (SNODE, see ishock_elm.h), including following message:
 * [ID] [x, y] [CHILD_ISHOCK_NODE_ID] [BND_1_ID] [BND_2_ID] ...
 *
 * In ishock edge, the shock type can be:
 * 1: point-point
 * 2: point-line
 * 4: line-line
 * 8: contact
 * with the message:
 * [ID] [TYPE] [FROM_NODE_ID] [TO_NODE_ID] [LEFT_BND_ID] [RIGHT_BND_ID]
 * TODO: what is Eta for ishock edge? How to input Eta when creating an ishock edge?
 * TODO: when arc fitting, we need add other types of ishock edge
 *
 * After finishing .isf loader, there are two members for accessing:
 * 1. _ishock_graph
 * 2. _belms
 * which means, there is no boundary info in _ishock_graph, only ishock nodes and edges
 * We create an boundary elemnts list for accessing the relative boundaries (i.e. _belms)
 */

#include <clocale>
#include "dbsk2d_load_isf_process.h"
#include "dbsk2d/dbsk2d_boundary.h"
#include "dbsk2d/dbsk2d_ishock_pointpoint.h"
#include "dbsk2d/dbsk2d_ishock_pointline.h"
#include "dbsk2d/dbsk2d_ishock_lineline.h"
#include "dbsk2d/dbsk2d_ishock_contact.h"

dbsk2d_load_isf_process::dbsk2d_load_isf_process() : bpro1_process()
{
    _ishock_graph = new dbsk2d_ishock_graph();
    if( !parameters()->add( "Input file <filename...>" , "-isfinput" , bpro1_filepath("ishock.isf","*.isf") ) )
    {
        vcl_cerr << "ERROR: Adding parameters in " __FILE__ << vcl_endl;
    }
}

bpro1_process *dbsk2d_load_isf_process::clone() const
{
    return new dbsk2d_load_isf_process(*this);
}

vcl_vector< vcl_string > dbsk2d_load_isf_process::get_input_type()
{
    vcl_vector< vcl_string > to_return;
    return to_return;
}

vcl_vector< vcl_string > dbsk2d_load_isf_process::get_output_type()
{
    vcl_vector< vcl_string > to_return;
    to_return.push_back( "shock" );
    return to_return;
}

bool dbsk2d_load_isf_process::execute()
{
    bpro1_filepath input;
    parameters()->get_value("-isfinput", input);
    vcl_string input_file = input.path;

    // set _ishock_graph & _belms in this function
    if(!loadISF(input_file))
    {
        vcl_cout << "ERROR: FAILED IN FUNC loadISF(input_file)." << vcl_endl;
        return false;
    }

    return true;
}

bool dbsk2d_load_isf_process::loadISF(vcl_string file_name)
{
    vcl_ifstream infp(file_name.c_str(), vcl_ios::in);
    int mode = 0;
    char start;

    // for bnd
    int bnd_id, bnd_type;
    double x1, y1, x2, y2;

    // for ishock node
    int node_id, child_ishock_node_id;
    double node_x, node_y;
    vcl_vector<int> bnd_ids;

    // for ishock edge
    int edge_id, edge_type, from_id, to_id, left_id, right_id;


    if (!infp){
        vcl_cerr << " Error opening file  " << file_name.c_str() << vcl_endl;
        return false;
    }

    char buf[512];
    int line_count = 0;

    while(!infp.eof())
    {
        infp.getline(buf, 512);

        if (buf[0] == '#')
        {
            mode = 0;
            continue;
        }
        else if (!vcl_strncmp(buf, "[BOUNDARY_BEGIN]", sizeof("[BOUNDARY_BEGIN]")-1))
        {
            mode = 1;
            continue;
        }
        else if (!vcl_strncmp(buf, "[ISHOCK_NODE_BEGIN]", sizeof("[ISHOCK_NODE_BEGIN]")-1))
        {
            mode = 2;
            continue;
        }
        else if (!vcl_strncmp(buf, "[ISHOCK_EDGE_BEGIN]", sizeof("[ISHOCK_EDGE_BEGIN]")-1))
        {
            vcl_cout << "change3" << vcl_endl;
            mode = 3;
            continue;
        }
        else;

        if (mode == 0)
        {
            continue;
        }
        else if (mode == 1)
        {
            vcl_cout << "read a boundary" << vcl_endl;
            // read a boundary
            sscanf(buf, "%d %d", &bnd_id, &bnd_type);
            if (bnd_type == 0)
            {
                // read bpoint
                sscanf(buf, "%d %d %lf %lf", &bnd_id, &bnd_type, &x1, &y1);
                dbsk2d_ishock_bpoint* new_bpoint = new dbsk2d_ishock_bpoint(x1, y1, bnd_id);

                _belms.push_back(new_bpoint);
            }
            else if (bnd_type == 1)
            {
                // read bline
                sscanf(buf, "%d %d %lf %lf %lf %lf", &bnd_id, &bnd_type, &x1, &y1, &x2, &y2);

                // create bline with p1 & p2
                dbsk2d_ishock_bpoint* start_bpoint = new dbsk2d_ishock_bpoint(x1, y1);
                dbsk2d_ishock_bpoint* end_bpoint = new dbsk2d_ishock_bpoint(x2, y2);
                dbsk2d_ishock_bline* new_bline = new dbsk2d_ishock_bline(start_bpoint, end_bpoint, bnd_id);

                _belms.push_back(new_bline);
            }
            else if (bnd_type == 2)
            {
                //todo barc
            }
            else
            {
                vcl_cout << "ERROR: NO SUCH BOUNDARY TYPE." << vcl_endl;
            }
        }
        else if (mode == 2)
        {
            vcl_cout << "read an ishock node" << vcl_endl;
            // read an ishock node
            sscanf(buf, "%d %lf %lf %d", &node_id, &node_x, &node_y, &child_ishock_node_id);
            get_bnds(buf, bnd_ids);

            // create a new ishock node
            vgl_point_2d<double> pt = vgl_point_2d<double>(node_x, node_y);
            dbsk2d_ishock_node* new_node = new dbsk2d_ishock_node(node_id, 0, pt);

            // use bnd ids to find out boundary that form the node
            for (int i = 0; i < bnd_ids.size(); i++)
            {
                new_node->add_belm(get_bnd_from_index(bnd_ids[i]));
                vcl_cout << "here" << vcl_endl;
            }

            _ishock_graph->add_vertex(new_node);
            bnd_ids.clear();
        }
        else if (mode == 3)
        {
            vcl_cout << "read an ishock edge" << vcl_endl;
            vcl_cout << buf << vcl_endl;
            // read an ishock edge
            sscanf(buf, "%d %d %d %d %d %d", &edge_id, &edge_type, &from_id, &to_id, &left_id, &right_id);

            // create a new ishock edge
            dbsk2d_ishock_node* p_shock_node = get_ishock_node_from_index(from_id);
            dbsk2d_ishock_node* c_shock_node = get_ishock_node_from_index(to_id);
            dbsk2d_ishock_belm* l_bnd_edge = get_bnd_from_index(left_id);
            dbsk2d_ishock_belm* r_bnd_edge = get_bnd_from_index(right_id);

            // this edge is point-point
            if (edge_type == 1)
            {
                dbsk2d_ishock_pointpoint* new_pp_edge = new dbsk2d_ishock_pointpoint(edge_id, 0, p_shock_node, l_bnd_edge, r_bnd_edge, 0, 0);

                // set relation between node and edge
                if (p_shock_node)
                    p_shock_node->set_cShock(new_pp_edge);
                new_pp_edge->set_cSNode(c_shock_node);

                _ishock_graph->add_edge(new_pp_edge);
            }
            // this edge is point-line
            else if (edge_type == 2)
            {
                dbsk2d_ishock_pointline* new_pl_edge = new dbsk2d_ishock_pointline(edge_id, 0, p_shock_node, l_bnd_edge, r_bnd_edge, 0, 0);

                // set relation between node and edge
                if (p_shock_node)
                    p_shock_node->set_cShock(new_pl_edge);
                new_pl_edge->set_cSNode(c_shock_node);

                _ishock_graph->add_edge(new_pl_edge);
            }
            // this edge is line-line
            else if (edge_type == 4)
            {
                dbsk2d_ishock_lineline* new_ll_edge = new dbsk2d_ishock_lineline(edge_id, 0, p_shock_node, l_bnd_edge, r_bnd_edge, 0, 0);

                // set relation between node and edge
                if (p_shock_node)
                    p_shock_node->set_cShock(new_ll_edge);
                new_ll_edge->set_cSNode(c_shock_node);

                _ishock_graph->add_edge(new_ll_edge);
            }
            // this edge is contact shock (only has a point)
            else if (edge_type == 8)
            {
                // todo: some contact is very confused, their from_node is -1 as well as to_node ( is it an unlimited line? )
                // we ignore such type of contact
                if (!c_shock_node)
                    continue;

                vgl_point_2d<double> pt = c_shock_node->origin();
                dbsk2d_ishock_contact* new_c_edge = new dbsk2d_ishock_contact(edge_id, l_bnd_edge, r_bnd_edge, pt, 0, 0, 0, 0, 0);

                // set relation between node and edge
                new_c_edge->set_cSNode(c_shock_node);

                _ishock_graph->add_edge(new_c_edge);
            }
            else
            {
                // todo : since we use line-fitting in current stage, we do not consider point-arc, arc-arc ...
                // currently, we only use point-point, point-line, line-line and contact type of ishock edge
                vcl_cout << "ERROR: NO SUCH TYPE OF ISHOCK EDGE." << vcl_endl;
            }
        }
        else
        {
            vcl_cout << "ERROR: NO SUCH MODE." << vcl_endl;
            return false;
        }
    }

    infp.close();
    return true;
}

void dbsk2d_load_isf_process::get_bnds(char *buf, vcl_vector<int> &bnd_ids)
{
    int num_bnds;
    int num_blanks = 0;
    int bnd_id1, bnd_id2, bnd_id3;
    float dummy;

    for (int i = 0; i < vcl_strlen(buf); i++)
    {
        if (buf[i] == ' ')
            num_blanks++;
    }
    // the num of bnds is not a const
    // pay attention to the blank at the end of each node line
    num_bnds = num_blanks - 4;

    // typically, there are 2 or 3 bnds for generating a ishock node
    if (num_bnds == 2)
    {
        sscanf(buf, "%d %f %f %d %d %d", &dummy, &dummy, &dummy, &dummy, &bnd_id1, &bnd_id2);
        bnd_ids.push_back(bnd_id1);
        bnd_ids.push_back(bnd_id2);
    }
    else if (num_bnds == 3)
    {
        sscanf(buf, "%d %f %f %d %d %d %d", &dummy, &dummy, &dummy, &dummy, &bnd_id1, &bnd_id2, &bnd_id3);
        bnd_ids.push_back(bnd_id1);
        bnd_ids.push_back(bnd_id2);
        bnd_ids.push_back(bnd_id3);
    }
    else
    {
        vcl_cout << buf << vcl_endl;
        vcl_cout << "ERROR: COUNT BNDS IN ISHOCK NODE (NEITHER 2 NOR 3)." << vcl_endl;
    }
}

dbsk2d_ishock_node *dbsk2d_load_isf_process::get_ishock_node_from_index(int index)
{
    if (index == -1)
        return NULL;
    else
    {
        ishock_node_list nodes = _ishock_graph->all_nodes();
        for (ishock_node_list_iter inl = nodes.begin(); inl != nodes.end() ; inl++)
        {
            if ((*inl)->id() == index)
                return (*inl);
        }
        vcl_cout << "ERROR: ISHOCK NODE NOT FOUND." << vcl_endl;
        return NULL;
    }
}

dbsk2d_ishock_belm *dbsk2d_load_isf_process::get_bnd_from_index(int index)
{
    if (index == -1)
        return NULL;
    else
    {
        for (vcl_vector< dbsk2d_ishock_belm* >::iterator it = _belms.begin(); it != _belms.end(); it++)
        {
            if ((*it)->id() == index)
                return (*it);
        }
        vcl_cout << "ERROR: BOUNDARY ELEMENT NOT FOUND. Index : " << index << vcl_endl;
        return NULL;
    }
}
