//
// Created by wenhan on 9/4/17.
//

#include "dbsk2d_load_osf_process.h"

#include <clocale>
#include "dbsk2d_load_osf_process.h"
#include "dbsk2d/dbsk2d_boundary.h"
#include "dbsk2d/dbsk2d_ishock_pointpoint.h"
#include "dbsk2d/dbsk2d_ishock_pointline.h"
#include "dbsk2d/dbsk2d_ishock_lineline.h"
#include "dbsk2d/dbsk2d_ishock_contact.h"
#include "dbsk2d/pro/dbsk2d_shock_storage.h"
#include "dbsk2d/algo/dbsk2d_prune_ishock.h"

dbsk2d_load_osf_process::dbsk2d_load_osf_process() : bpro1_process()
{
    _ishock_graph = new dbsk2d_ishock_graph();
    _output_shock = dbsk2d_shock_storage_new();
    _output_shock->set_name("wenhanshi");
    if( !parameters()->add( "Input file <filename...>" , "-osfinput" , bpro1_filepath("output_shock.osf","*.osf") )
        || !parameters()->add( "Prune Threshold" , "-threshold" , (float)1.0 ))
    {
        vcl_cerr << "ERROR: Adding parameters in " __FILE__ << vcl_endl;
    }
}

bpro1_process *dbsk2d_load_osf_process::clone() const
{
    return new dbsk2d_load_osf_process(*this);
}

vcl_vector< vcl_string > dbsk2d_load_osf_process::get_input_type()
{
    vcl_vector< vcl_string > to_return;
    return to_return;
}

vcl_vector< vcl_string > dbsk2d_load_osf_process::get_output_type()
{
    vcl_vector< vcl_string > to_return;
    to_return.push_back( " " );
    return to_return;
}

bool dbsk2d_load_osf_process::execute()
{
    bpro1_filepath input;
    parameters()->get_value("-osfinput", input);
    vcl_string input_file = input.path;

    // set _ishock_graph & _belms in this function
    if(!load_osf(input_file))
    {
        vcl_cout << "ERROR: FAILED IN FUNC load_osf(input_file)." << vcl_endl;
        return false;
    }
    // let belms be in the boundar
    dbsk2d_boundary_sptr boundary = new dbsk2d_boundary();
    boundary->set_belms(_belms);

    // very important, to compute in information
    _ishock_graph->update_shocks();

    // generate coarse ishock by pruner
    dbsk2d_shock_graph_sptr coarse_ishock = new dbsk2d_shock_graph();
    float prune_threshold = 1.0;
    dbsk2d_prune_ishock ishock_pruner(_ishock_graph, coarse_ishock);
    parameters()->get_value( "-threshold" , prune_threshold );
    ishock_pruner.prune(prune_threshold);
    ishock_pruner.compile_coarse_shock_graph();

    // as the original design, output_shock should include three useful items:
    // 1. boundary (only belms)
    // 2. ishock graph (nodes and links)
    // 3. coarse shock graph (formed by pruner when loading ishock graph)
    _output_shock->set_boundary(boundary);
    _output_shock->set_ishock_graph(_ishock_graph);
    _output_shock->set_shock_graph(coarse_ishock);

    return true;
}

bool dbsk2d_load_osf_process::load_osf(vcl_string file_name)
{
    vcl_ifstream infp(file_name.c_str(), vcl_ios::in);
    int mode = 0;
    char start;
    double dummy;

    // for bnd
    int bnd_id, bnd_type, start_id, end_id, num_bnd, bnd_id1, bnd_id2, bnd_id3, bnd_id4;
    double x1, y1, x2, y2, u, n, l;

    // for ishock node
    int node_id, child_shock1, child_shock2, num_pshock, pshock_1, pshock_2, pshock_3, pshock_4;
    double node_x, node_y, start_time, end_time;
    vcl_vector<int> bnd_ids;

    // for ishock edge
    int edge_id, edge_type, from_id, to_id, left_id, right_id, nu;
    double ls_eta, rs_eta, ls_tau, le_tau, rs_tau, re_tau, h, edge_n, edge_u, ur, phi, l_delta, r_delta, ul, sigma, theta_l, theta_r, ll, lr;


    if (!infp){
        vcl_cerr << " Error opening file  " << file_name.c_str() << vcl_endl;
        return false;
    }

    char buf[512];
    int line_count = 0;

    while(!infp.eof())
    {
        infp.getline(buf, 512);

        // the end of the line
        if (buf[0] == '\0')
            continue;

        if (buf[0] == '#')
        {
            mode = 0;
            continue;
        }
        else if (!vcl_strncmp(buf, "[BOUNDARY_POINT_BEGIN]", sizeof("[BOUNDARY_POINT_BEGIN]")-1))
        {
            mode = 1;
            continue;
        }
        else if (!vcl_strncmp(buf, "[BOUNDARY_LINE_BEGIN]", sizeof("[BOUNDARY_LINE_BEGIN]")-1))
        {
            mode = 2;
            continue;
        }
        else if (!vcl_strncmp(buf, "[ISHOCK_NODE_BEGIN]", sizeof("[ISHOCK_NODE_BEGIN]")-1))
        {
            mode = 3;
            continue;
        }
        else if (!vcl_strncmp(buf, "[ISHOCK_EDGE_BEGIN]", sizeof("[ISHOCK_EDGE_BEGIN]")-1))
        {
            mode = 4;
            continue;
        }
        else if(!vcl_strncmp(buf, "[OTHER_INFO_BEGIN]", sizeof("[OTHER_INFO_BEGIN]")-1))
        {
            mode = 5;
            continue;
        }
        else;

        if (mode == 0)
        {
            continue;
        }
        else if (mode == 1)
        {
            //vcl_cout << "read a boundary point" << vcl_endl;
            // read a bpoint
            sscanf(buf, "%d %d %lf %lf",
                   &bnd_id,
                   &bnd_type,
                   &x1,
                   &y1
            );
            dbsk2d_ishock_bpoint* new_bpoint = new dbsk2d_ishock_bpoint(x1, y1, bnd_id);
            _belms.push_back(new_bpoint);
        }
        else if (mode == 2)
        {
            //vcl_cout << "read a boundary line" << vcl_endl;
            // read a bline
            sscanf(buf, "%d %d %d %d %lf %lf %lf",
                   &bnd_id,
                   &bnd_type,
                   &start_id,
                   &end_id,
                   &u,
                   &n,
                   &l
            );
            dbsk2d_ishock_bpoint* start_bpoint = (dbsk2d_ishock_bpoint*)get_bnd_from_index(start_id);
            dbsk2d_ishock_bpoint* end_bpoint = (dbsk2d_ishock_bpoint*)get_bnd_from_index(end_id);
            // create bline with p1 & p2
            dbsk2d_ishock_bline* new_bline = new dbsk2d_ishock_bline(start_bpoint, end_bpoint, bnd_id);
            new_bline->set_u(u);
            new_bline->set_n(n);
            new_bline->set_l(l);
            _belms.push_back(new_bline);
        }
        else if (mode == 3)
        {
            //vcl_cout << "read an ishock node" << vcl_endl;
            // read an ishock node
            sscanf(buf, "%d %lf %lf %lf %lf %d %d %d",
                   &node_id,
                   &node_x,
                   &node_y,
                   &start_time,
                   &end_time,
                   &child_shock1,
                   &child_shock2,
                   &num_bnd
            );
            if (num_bnd == 2)
            {
                sscanf(buf, "%d %f %f %f %f %d %d %d %d %d",
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &bnd_id1,
                       &bnd_id2
                );
                bnd_ids.push_back(bnd_id1);
                bnd_ids.push_back(bnd_id2);
            }
            else if (num_bnd == 3)
            {
                sscanf(buf, "%d %f %f %f %f %d %d %d %d %d %d",
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &bnd_id1,
                       &bnd_id2,
                       &bnd_id3
                );
                bnd_ids.push_back(bnd_id1);
                bnd_ids.push_back(bnd_id2);
                bnd_ids.push_back(bnd_id3);
            }
            else if (num_bnd == 4)
            {
                sscanf(buf, "%d %f %f %f %f %d %d %d %d %d %d %d",
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &dummy,
                       &bnd_id1,
                       &bnd_id2,
                       &bnd_id3,
                       &bnd_id4
                );
                bnd_ids.push_back(bnd_id1);
                bnd_ids.push_back(bnd_id2);
                bnd_ids.push_back(bnd_id3);
                bnd_ids.push_back(bnd_id4);
            }
            else
            {
                vcl_cout << buf << vcl_endl;
                vcl_cout << "ERROR: COUNT BNDS IN ISHOCK NODE (!= 2 or 3 or 4)." << vcl_endl;
            }

            // create a new ishock node
            vgl_point_2d<double> pt = vgl_point_2d<double>(node_x, node_y);
            dbsk2d_ishock_node* new_node = new dbsk2d_ishock_node(node_id, start_time, pt);
            new_node->setEndTime(end_time);
            new_node->set_cshock_id_1(child_shock1);
            new_node->set_cshock_id_2(child_shock2);
            new_node->ex_pts().push_back(pt);

            // use bnd ids to find out boundary that form the node
            for (int i = 0; i < bnd_ids.size(); i++)
            {
                new_node->add_belm(get_bnd_from_index(bnd_ids[i]));
            }

            _ishock_graph->add_vertex(new_node);
            bnd_ids.clear();
        }
        else if (mode == 4)
        {
            //vcl_cout << "read an ishock edge" << vcl_endl;
            //vcl_cout << buf << vcl_endl;
            // read an ishock edge
            sscanf(buf, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d",
                   &edge_id,
                   &edge_type,
                   &start_time,
                   &end_time,
                   &ls_eta,
                   &rs_eta,
                   &ls_tau,
                   &le_tau,
                   &rs_tau,
                   &re_tau,
                   &h,
                   &from_id,
                   &to_id,
                   &left_id,
                   &right_id
            );

            // create a new ishock edge
            dbsk2d_ishock_node* p_shock_node = get_ishock_node_from_index(from_id);
            dbsk2d_ishock_node* c_shock_node = get_ishock_node_from_index(to_id);
            dbsk2d_ishock_belm* l_bnd_edge = get_bnd_from_index(left_id);
            dbsk2d_ishock_belm* r_bnd_edge = get_bnd_from_index(right_id);

            // this edge is point-point
            if (edge_type == 1)
            {
                sscanf(buf, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %lf %lf",
                       &edge_id,
                       &edge_type,
                       &start_time,
                       &end_time,
                       &ls_eta,
                       &rs_eta,
                       &ls_tau,
                       &le_tau,
                       &rs_tau,
                       &re_tau,
                       &h,
                       &from_id,
                       &to_id,
                       &left_id,
                       &right_id,
                       &edge_n,
                       &edge_u
                );
                dbsk2d_ishock_pointpoint* new_pp_edge = new dbsk2d_ishock_pointpoint(edge_id,
                                                                                     start_time,
                                                                                     p_shock_node,
                                                                                     l_bnd_edge,
                                                                                     r_bnd_edge,
                                                                                     ls_eta,
                                                                                     rs_eta);

                // set relation between node and edge
                if (p_shock_node)
                    p_shock_node->set_cShock(new_pp_edge);
                new_pp_edge->set_cSNode(c_shock_node);
                new_pp_edge->set_lbnd(l_bnd_edge);
                new_pp_edge->set_rbnd(r_bnd_edge);
                new_pp_edge->setEndTime(end_time);
                new_pp_edge->setLsTau(ls_tau);
                new_pp_edge->setLeTau(le_tau);
                new_pp_edge->setRsTau(rs_tau);
                new_pp_edge->setReTau(re_tau);
                new_pp_edge->set_H(h);
                new_pp_edge->set_n(edge_n);
                new_pp_edge->set_u(edge_u);

                _ishock_graph->add_edge(new_pp_edge);
            }
                // this edge is point-line
            else if (edge_type == 2)
            {
                // need to re-read for _nu
                sscanf(buf, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %lf %lf %lf %lf %lf",
                       &edge_id,
                       &edge_type,
                       &start_time,
                       &end_time,
                       &ls_eta,
                       &rs_eta,
                       &ls_tau,
                       &le_tau,
                       &rs_tau,
                       &re_tau,
                       &h,
                       &from_id,
                       &to_id,
                       &left_id,
                       &right_id,
                       &nu,
                       &edge_u,
                       &edge_n,
                       &l_delta,
                       &r_delta,
                       &l
                );

                dbsk2d_ishock_pointline* new_pl_edge = new dbsk2d_ishock_pointline(edge_id,
                                                                                   start_time,
                                                                                   p_shock_node,
                                                                                   l_bnd_edge,
                                                                                   r_bnd_edge,
                                                                                   ls_eta,
                                                                                   rs_eta);

                // set relation between node and edge
                if (p_shock_node)
                    p_shock_node->set_cShock(new_pl_edge);
                new_pl_edge->set_cSNode(c_shock_node);
                new_pl_edge->set_lbnd(l_bnd_edge);
                new_pl_edge->set_rbnd(r_bnd_edge);
                new_pl_edge->setEndTime(end_time);
                new_pl_edge->setLsTau(ls_tau);
                new_pl_edge->setLeTau(le_tau);
                new_pl_edge->setRsTau(rs_tau);
                new_pl_edge->setReTau(re_tau);
                new_pl_edge->set_H(h);
                new_pl_edge->set_nu(nu);
                new_pl_edge->set_u(edge_u);
                new_pl_edge->set_n(edge_n);
                new_pl_edge->set_ldelta(l_delta);
                new_pl_edge->set_rdelta(r_delta);
                new_pl_edge->set_l(l);

                _ishock_graph->add_edge(new_pl_edge);
            }
                // this edge is line-line
            else if (edge_type == 4)
            {
                sscanf(buf, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf",
                       &edge_id,
                       &edge_type,
                       &start_time,
                       &end_time,
                       &ls_eta,
                       &rs_eta,
                       &ls_tau,
                       &le_tau,
                       &rs_tau,
                       &re_tau,
                       &h,
                       &from_id,
                       &to_id,
                       &left_id,
                       &right_id,
                       &ur,
                       &phi,
                       &ul,
                       &sigma,
                       &theta_l,
                       &theta_r,
                       &ll,
                       &lr
                );
                dbsk2d_ishock_lineline* new_ll_edge = new dbsk2d_ishock_lineline(edge_id,
                                                                                 start_time,
                                                                                 p_shock_node,
                                                                                 l_bnd_edge,
                                                                                 r_bnd_edge,
                                                                                 ls_eta,
                                                                                 rs_eta);

                // set relation between node and edge
                if (p_shock_node)
                    p_shock_node->set_cShock(new_ll_edge);
                new_ll_edge->set_cSNode(c_shock_node);
                new_ll_edge->set_lbnd(l_bnd_edge);
                new_ll_edge->set_rbnd(r_bnd_edge);
                new_ll_edge->setEndTime(end_time);
                new_ll_edge->setLsTau(ls_tau);
                new_ll_edge->setLeTau(le_tau);
                new_ll_edge->setRsTau(rs_tau);
                new_ll_edge->setReTau(re_tau);
                new_ll_edge->set_H(h);
                new_ll_edge->set_ur(ur);
                new_ll_edge->set_phi(phi);
                new_ll_edge->set_ul(ul);
                new_ll_edge->set_sigma(sigma);
                new_ll_edge->set_thetaL(theta_l);
                new_ll_edge->set_thetaR(theta_r);
                new_ll_edge->set_lL(ll);
                new_ll_edge->set_lR(lr);

                _ishock_graph->add_edge(new_ll_edge);
            }
                // this edge is contact shock (only has a point)
            else if (edge_type == 8)
            {
                sscanf(buf, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %lf",
                       &edge_id,
                       &edge_type,
                       &start_time,
                       &end_time,
                       &ls_eta,
                       &rs_eta,
                       &ls_tau,
                       &le_tau,
                       &rs_tau,
                       &re_tau,
                       &h,
                       &from_id,
                       &to_id,
                       &left_id,
                       &right_id,
                       &edge_n
                );
                // we ignore such type of contact
                if (!c_shock_node)
                    continue;
                vgl_point_2d<double> pt = c_shock_node->origin();
                dbsk2d_ishock_contact* new_c_edge = new dbsk2d_ishock_contact(edge_id,
                                                                              l_bnd_edge,
                                                                              r_bnd_edge,
                                                                              pt,
                                                                              0,
                                                                              0,
                                                                              0,
                                                                              ls_eta,
                                                                              rs_eta);

                // set relation between node and edge
                new_c_edge->set_cSNode(c_shock_node);
                new_c_edge->set_lbnd(l_bnd_edge);
                new_c_edge->set_rbnd(r_bnd_edge);
                new_c_edge->setEndTime(end_time);
                new_c_edge->setLsTau(ls_tau);
                new_c_edge->setLeTau(le_tau);
                new_c_edge->setRsTau(rs_tau);
                new_c_edge->setReTau(re_tau);
                new_c_edge->set_H(h);
                new_c_edge->set_n(edge_n);

                _ishock_graph->add_edge(new_c_edge);
            }
            else
            {
                // todo : since we use line-fitting in current stage, we do not consider point-arc, arc-arc ...
                // currently, we only use point-point, point-line, line-line and contact type of ishock edge
                vcl_cout << "ERROR: NO SUCH TYPE OF ISHOCK EDGE." << vcl_endl;
            }
        }
        else if (mode == 5)
        {
            sscanf(buf, "%d %d", &node_id, &num_pshock);
            if (num_pshock == 0);
            else if (num_pshock == 1)
            {
                dbsk2d_ishock_node* node = get_ishock_node_from_index(node_id);
                sscanf(buf, "%d %d %d", &node_id, &num_pshock, &pshock_1);
                node->add_pShock(get_ishock_edge_from_index(pshock_1));
            }
            else if (num_pshock == 2)
            {
                dbsk2d_ishock_node* node = get_ishock_node_from_index(node_id);
                sscanf(buf, "%d %d %d %d", &node_id, &num_pshock, &pshock_1, &pshock_2);
                node->add_pShock(get_ishock_edge_from_index(pshock_1));
                node->add_pShock(get_ishock_edge_from_index(pshock_2));
            }
            else if (num_pshock == 3)
            {
                dbsk2d_ishock_node* node = get_ishock_node_from_index(node_id);
                sscanf(buf, "%d %d %d %d %d", &node_id, &num_pshock, &pshock_1, &pshock_2, &pshock_3);
                node->add_pShock(get_ishock_edge_from_index(pshock_1));
                node->add_pShock(get_ishock_edge_from_index(pshock_2));
                node->add_pShock(get_ishock_edge_from_index(pshock_3));
            }
            else if (num_pshock == 4)
            {
                dbsk2d_ishock_node* node = get_ishock_node_from_index(node_id);
                sscanf(buf, "%d %d %d %d %d %d", &node_id, &num_pshock, &pshock_1, &pshock_2, &pshock_3, &pshock_4);
                node->add_pShock(get_ishock_edge_from_index(pshock_1));
                node->add_pShock(get_ishock_edge_from_index(pshock_2));
                node->add_pShock(get_ishock_edge_from_index(pshock_3));
                node->add_pShock(get_ishock_edge_from_index(pshock_4));
            }
            else
            {
                vcl_cout << buf << vcl_endl;
                vcl_cout << "ERROR: IMPOSSIBLE NUM OF PSHOCKS (!= 2, 3 or 4)." << vcl_endl;
            }
        }
        else
        {
            vcl_cout << "ERROR: NO SUCH MODE." << vcl_endl;
            return false;
        }
    }

    infp.close();

    // add all child shock links (# of links is 1 or 2) for each ishock node
    for (ishock_node_list_iter inl = _ishock_graph->all_nodes().begin(); inl != _ishock_graph->all_nodes().end() ; inl++)
    {
        (*inl)->set_cShock(get_ishock_edge_from_index((*inl)->get_cshock_id_1()));
        (*inl)->set_cShock2(get_ishock_edge_from_index((*inl)->get_cshock_id_2()));
    }

    return true;
}

dbsk2d_ishock_node *dbsk2d_load_osf_process::get_ishock_node_from_index(int index)
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

dbsk2d_ishock_belm *dbsk2d_load_osf_process::get_bnd_from_index(int index)
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

dbsk2d_ishock_edge *dbsk2d_load_osf_process::get_ishock_edge_from_index(int index)
{
    if (index == -1)
        return NULL;
    else
    {
        ishock_edge_list edges = _ishock_graph->all_edges();
        for (ishock_edge_list_iter iel = edges.begin(); iel != edges.end() ; iel++)
        {
            if ((*iel)->id() == index)
                return (*iel);
        }
        vcl_cout << "ERROR: ISHOCK NODE NOT FOUND." << vcl_endl;
        return NULL;
    }
}
