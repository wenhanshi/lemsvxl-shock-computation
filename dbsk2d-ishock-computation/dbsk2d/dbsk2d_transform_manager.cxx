// This is brcv/shp/dbsk2d/algo/dbsk2d_transform_manager.cxx

//:
// \file
#include "dbsk2d_transform_manager.h"
#include <vil3d/vil3d_load.h>
#include <vil/vil_load.h>
#include <vil/vil_bilin_interp.h>
#include <vil3d/vil3d_tricub_interp.h>
#include <vcl_fstream.h>

#include <vil/vil_plane.h>
#include <vil/vil_save.h>
#include <vil/vil_new.h>
#include <vnl/algo/vnl_chi_squared.h>
#include <vnl/algo/vnl_svd.h>

#include <vil/vil_convert.h>
#include <vsol/vsol_box_2d.h>
#include <vsol/vsol_polygon_2d.h>
#include <vsol/vsol_polygon_2d_sptr.h>

#include <vcl_cmath.h>
#include <bsol/bsol_algs.h>

#include <vgl/vgl_polygon_scan_iterator.h>
#include <vgl/algo/vgl_convex_hull_2d.h>

#include "dbsk2d_bnd_contour.h"

#include <vgl/algo/vgl_fit_lines_2d.h>
#include <vsol/vsol_polyline_2d.h>

#include "../dbsol/dbsol_curve_algs.h"

#include "dbsk2d_shock_graph_sptr.h"
#include "dbsk2d_xshock_edge.h"
#include "algo/dbsk2d_sample_ishock.h"
#include "algo/dbsk2d_ishock_grouping_transform.h"
#include <bbas/bsta/bsta_histogram.h>

#include <vgl/vgl_intersection.h>
#include <vgl/vgl_line_2d.h>

#include "../vidpro1/storage/vidpro1_vsol2D_storage_sptr.h"
#include "../vidpro1/storage/vidpro1_vsol2D_storage.h"
#include <vsol/vsol_line_2d.h>

#include "../dbsk2d_compute_ishock_process.h"
#include "pro/dbsk2d_shock_storage.h"
#include <vgl/vgl_area.h>

dbsk2d_transform_manager::dbsk2d_transform_manager()
    :image_(0),
     threshold_(0),
     dist_volume_(0),
     gPb_image_(0),
     out_folder_(),
     out_prefix_(),
     logistic_beta0_(1.0),
     logistic_beta1_(0.0),
     id_(0),
     normalization_(1.0),
     diag_(1.0)

{


}

dbsk2d_transform_manager::~dbsk2d_transform_manager()
{
    destroy_singleton();
}

void dbsk2d_transform_manager::destroy_singleton()
{
   
    image_=0;
}

void dbsk2d_transform_manager::read_in_training_data(vcl_string filename)
{

    vcl_ifstream file (filename.c_str(), 
                       vcl_ios::in|vcl_ios::binary|vcl_ios::ate);
    double* memblock(0);
    if (file.is_open())
    {
        vcl_ifstream::pos_type size = file.tellg();
        memblock = new double[size/sizeof(double)];
        file.seekg (0, vcl_ios::beg);
        file.read ((char *) memblock, size);
        file.close();

        // Read in dimensions
        unsigned int ni=memblock[0];
        unsigned int nj=memblock[1];
        unsigned int nk=memblock[2];

        vcl_cout<<"Reading in a "<<ni<<" by "<<nj<<" by "<< nk 
                <<" volume of distances"
                <<vcl_endl;

        dist_volume_.set_size(ni,nj,nk,1);
        unsigned int index=3;

        for (unsigned k=0;k<nk;++k)
        {
            for (unsigned j=0;j<nj;++j)
            {
                for (unsigned i=0;i<ni;++i)
                {
                    double value = memblock[index];
                    dist_volume_(i,j,k,0)=value;
                    index++;
                }
            }
        }

        delete[] memblock;
        memblock=0;
    }

    
}

void dbsk2d_transform_manager::read_in_gpb_data(vcl_string filename)
{

    vcl_ifstream file (filename.c_str(), 
                       vcl_ios::in|vcl_ios::binary|vcl_ios::ate);
    double* memblock(0);
    double max_gPb_value=0.0;
    if (file.is_open())
    {
        vcl_ifstream::pos_type size = file.tellg();
        memblock = new double[size/sizeof(double)];
        file.seekg (0, vcl_ios::beg);
        file.read ((char *) memblock, size);
        file.close();

        // Read in dimensions
        unsigned int ni=memblock[0];
        unsigned int nj=memblock[1];

        vcl_cout<<"Reading in a "<<ni<<" by "<<nj<<" gPb values"
                <<vcl_endl;

        gPb_image_.set_size(ni,nj);
        unsigned int index=2;

        for (unsigned j=0;j<nj;++j)
        {
            for (unsigned i=0;i<ni;++i)
            {
                double value = memblock[index];
                gPb_image_(i,j)=value;
                if ( value > max_gPb_value )
                {
                    max_gPb_value=value;
                }
                index++;
            }
        }
        

        delete[] memblock;
        memblock=0;
    }

    normalization_=max_gPb_value;
}

void dbsk2d_transform_manager::read_in_texton_data(vcl_string filename)
{

    vcl_ifstream file (filename.c_str(), 
                       vcl_ios::in|vcl_ios::binary|vcl_ios::ate);
    double* memblock(0);

    if (file.is_open())
    {
        vcl_ifstream::pos_type size = file.tellg();
        memblock = new double[size/sizeof(double)];
        file.seekg (0, vcl_ios::beg);
        file.read ((char *) memblock, size);
        file.close();

        // Read in dimensions
        unsigned int ni=memblock[0];
        unsigned int nj=memblock[1];

        vcl_cout<<"Reading in a "<<ni<<" by "<<nj<<" texton values"
                <<vcl_endl;

        texton_image_.set_size(ni,nj);
        unsigned int index=2;

        for (unsigned j=0;j<nj;++j)
        {
            for (unsigned i=0;i<ni;++i)
            {
                double value = memblock[index];
                texton_image_(i,j)=value+1.0;
                index++;
            }
        }
        

        delete[] memblock;
        memblock=0;
    }

}

double dbsk2d_transform_manager::contour_gpb_value(
    vcl_vector<dbsk2d_ishock_belm*>& belms)
{
    vcl_map<int, vgl_point_2d<double> > output_points;
  
    double perimeter=0.0;
    vcl_vector<dbsk2d_ishock_belm*>::iterator lit;  
    for (lit = belms.begin() ; lit != belms.end() ; ++lit)
    {
        dbsk2d_ishock_bline* bline = (dbsk2d_ishock_bline*)(*lit);
        double distance=vgl_distance(bline->s_pt()->pt(),
                                     bline->e_pt()->pt());

        output_points[bline->s_pt()->id()]=bline->s_pt()->pt();
        output_points[bline->e_pt()->id()]=bline->e_pt()->pt();

        perimeter=perimeter+distance;
    }     


    double summation=0.0;
    vcl_map<int,vgl_point_2d<double> >::iterator it;
    for ( it = output_points.begin() ; it != output_points.end(); ++it)
    {
        double x=(*it).second.x();
        double y=(*it).second.y();

        double gPb = vil_bilin_interp_safe_extend(gPb_image_,
                                                  x,
                                                  y);

        summation = summation + gPb;
    }
    


    return summation/perimeter;
}


double dbsk2d_transform_manager::region_gpb_value(
    vcl_vector<vgl_point_2d<double> >& grid)
{

    double summation=0.0;
    for ( unsigned int g=0; g < grid.size() ; ++g)
    {
        double x=grid[g].x();
        double y=grid[g].y();

        double gPb = vil_bilin_interp_safe_extend(gPb_image_,
                                                  x,
                                                  y);

        summation = summation + gPb;
    }

    return summation;
}

void dbsk2d_transform_manager::write_stats_closed(
    vgl_polygon<double>& polygon)
{


    // // Find one sheeted polygon
    // unsigned int f_index=0;
    // vgl_polygon<double> start_poly(polygon[f_index]);
    // double area=vgl_area(start_poly);

    // for (unsigned int s = 1; s < polygon.num_sheets(); ++s)
    // { 
        
    //     vgl_polygon<double> tempy(polygon[s]);
    //     double area_temp = vgl_area(tempy);
    //     if ( area_temp > area )
    //     {
    //          area = area_temp;
    //          f_index=s;

    //     }
        
    // }

    // vgl_polygon<double> model_poly(polygon[f_index]);

    // vidpro1_vsol2D_storage_sptr input_vsol = 
    //     vidpro1_vsol2D_storage_new();
    // vsol_polygon_2d_sptr region_poly = bsol_algs::poly_from_vgl(model_poly);
    // vsol_box_2d_sptr bbox=new vsol_box_2d();

    // // Enlarge bounding box from size
    // // Calculate xcenter, ycenter
    // double xcenter = image_->ni()/2.0;
    // double ycenter = image_->nj()/2.0;
    
    // // Translate to center and scale
    // double xmin_scaled = ((0-xcenter)*2)+xcenter;
    // double ymin_scaled = ((0-ycenter)*2)+ycenter;
    // double xmax_scaled = ((image_->ni()-xcenter)*2)+xcenter;
    // double ymax_scaled = ((image_->nj()-ycenter)*2)+ycenter;
    
    // bbox->add_point(xmin_scaled,ymin_scaled);
    // bbox->add_point(xmax_scaled,ymax_scaled);
        
    // vsol_polygon_2d_sptr box_poly = bsol_algs::poly_from_box(bbox);

    // input_vsol->add_object(box_poly->cast_to_spatial_object());
    // input_vsol->add_object(region_poly->cast_to_spatial_object());

    // /*********************** Shock Compute **********************************/
    // // Grab output from shock computation
    // vcl_vector<bpro1_storage_sptr> shock_results;
    // bool status=false;
    // {
    //     // 3) Create shock pro process and assign inputs 
    //     dbsk2d_compute_ishock_process shock_pro;

    //     shock_pro.clear_input();
    //     shock_pro.clear_output();
        
    //     shock_pro.add_input(0);
    //     shock_pro.add_input(input_vsol);

    //     // Set params
    //     status = shock_pro.execute();
    //     shock_pro.finish();

    //     // If ishock status is bad we will keep iterating with noise 
    //     // till we get a valid shock computation otherwise call it quits
    //     if (!status)
    //     {
    //         // Add noise to parameter set
    //         shock_pro.parameters()->set_value("-b_noise",true);
            
    //         // Clean up before we start running
    //         shock_pro.clear_input();
    //         shock_pro.clear_output();
            
    //         unsigned int i(0);
    //         unsigned int num_iterations = 5;
            
    //         for ( ; i < num_iterations; ++i)
    //         {
    //             vcl_cout<<vcl_endl;
    //             vcl_cout<<"************ Retry Compute Shock,iter: "
    //                     <<i+1<<" *************"<<vcl_endl;
                
    //             // Add inputs
    //             shock_pro.add_input(0);
    //             shock_pro.add_input(input_vsol);
                
    //             // Kick off process again
    //             status = shock_pro.execute();
    //             shock_pro.finish();
                
    //             if ( status )
    //             {
    //                 // We have produced valid shocks lets quit
    //                 break;
                    
    //             }
                
    //             // Clean up after ourselves
    //             shock_pro.clear_input();
    //             shock_pro.clear_output();
                
    //         }
    //     }

    //     if ( status )
    //     {
    //         shock_results = shock_pro.get_output();

    //         // Clean up after ourselves
    //         shock_pro.clear_input();
    //         shock_pro.clear_output();
            
    //     }
        
  
    // }

    // dbsk2d_shock_storage_sptr shock_storage;
    // shock_storage.vertical_cast(shock_results[0]);
  
    // dbsk2d_ishock_grouping_transform grouper(shock_storage
    //                                          ->get_ishock_graph());
    // grouper.grow_regions();


    // vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_node*> >
    //     fragments = grouper.get_outer_shock_nodes();
    // vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_edge*> >
    //     frag_edges = grouper.get_region_nodes();
    // vcl_map<unsigned int, vcl_vector<dbsk2d_ishock_belm*> > 
    //     frag_belms = grouper.get_region_belms();

    // vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_edge*> >::iterator it;
    // for ( it = frag_edges.begin() ; it != frag_edges.end() ; ++it)
    // {

    //     bool closed_region=(fragments[(*it).first].size()==0)?
    //         true:false;

    //     if ( closed_region && grouper.region_within_image((*it).first,-1))
    //     {
    //         break;
    //     }
    // }

    // vcl_vector<double> region_stats;
    // vcl_vector<double> app_stats;

    // // get polygon stats
    // grouper.get_region_stats((*it).first,
    //                          model_poly,region_stats);

    // this->get_appearance_stats(frag_edges[(*it).first],
    //                            frag_belms[(*it).first],
    //                            vgl_area(model_poly),
    //                            app_stats);

    // vcl_vector<double> total_stats;
    // total_stats.push_back(0.0); //depth , look at this again
    // total_stats.push_back(1.0); //path prob
    // total_stats.push_back(1.0); //region gap cost
    
    // for ( unsigned int p=0; p < region_stats.size() ; ++p)
    // {
    //     total_stats.push_back(region_stats[p]);
    // }
    
    // for ( unsigned int a=0; a < app_stats.size() ; ++a)
    // {
    //     total_stats.push_back(app_stats[a]);
    // }

    // this->write_output_region_stats(total_stats);

}


void dbsk2d_transform_manager::write_stats_closed(
    vcl_vector<dbsk2d_ishock_belm*>& belms)
{

    // vidpro1_vsol2D_storage_sptr input_vsol = 
    //     vidpro1_vsol2D_storage_new();

    // for ( unsigned int i=0; i < belms.size() ; ++i)
    // {
    //     if ( belms[i]->is_a_line() )
    //     {
    //         dbsk2d_ishock_bline* line_element = 
    //             dynamic_cast<dbsk2d_ishock_bline*>(belms[i]);

    //         // Add in contours for front 
    //         vsol_spatial_object_2d_sptr obj=
    //             new vsol_line_2d(line_element->s_pt()->pt(),
    //                              line_element->e_pt()->pt());
    //         input_vsol->add_object(obj);
    //     }
    // }

    // vsol_box_2d_sptr bbox=new vsol_box_2d();

    // // Enlarge bounding box from size
    // // Calculate xcenter, ycenter
    // double xcenter = image_->ni()/2.0;
    // double ycenter = image_->nj()/2.0;
    
    // // Translate to center and scale
    // double xmin_scaled = ((0-xcenter)*2)+xcenter;
    // double ymin_scaled = ((0-ycenter)*2)+ycenter;
    // double xmax_scaled = ((image_->ni()-xcenter)*2)+xcenter;
    // double ymax_scaled = ((image_->nj()-ycenter)*2)+ycenter;
    
    // bbox->add_point(xmin_scaled,ymin_scaled);
    // bbox->add_point(xmax_scaled,ymax_scaled);
        
    // vsol_polygon_2d_sptr box_poly = bsol_algs::poly_from_box(bbox);
    // input_vsol->add_object(box_poly->cast_to_spatial_object());

    // /*********************** Shock Compute **********************************/
    // // Grab output from shock computation
    // vcl_vector<bpro1_storage_sptr> shock_results;
    // bool status=false;
    // {
    //     // 3) Create shock pro process and assign inputs 
    //     dbsk2d_compute_ishock_process shock_pro;

    //     shock_pro.clear_input();
    //     shock_pro.clear_output();
        
    //     shock_pro.add_input(0);
    //     shock_pro.add_input(input_vsol);

    //     // Set params
    //     status = shock_pro.execute();
    //     shock_pro.finish();

    //     // If ishock status is bad we will keep iterating with noise 
    //     // till we get a valid shock computation otherwise call it quits
    //     if (!status)
    //     {
    //         // Add noise to parameter set
    //         shock_pro.parameters()->set_value("-b_noise",true);
            
    //         // Clean up before we start running
    //         shock_pro.clear_input();
    //         shock_pro.clear_output();
            
    //         unsigned int i(0);
    //         unsigned int num_iterations = 5;
            
    //         for ( ; i < num_iterations; ++i)
    //         {
    //             vcl_cout<<vcl_endl;
    //             vcl_cout<<"************ Retry Compute Shock,iter: "
    //                     <<i+1<<" *************"<<vcl_endl;
                
    //             // Add inputs
    //             shock_pro.add_input(0);
    //             shock_pro.add_input(input_vsol);
                
    //             // Kick off process again
    //             status = shock_pro.execute();
    //             shock_pro.finish();
                
    //             if ( status )
    //             {
    //                 // We have produced valid shocks lets quit
    //                 break;
                    
    //             }
                
    //             // Clean up after ourselves
    //             shock_pro.clear_input();
    //             shock_pro.clear_output();
                
    //         }
    //     }

    //     if ( status )
    //     {
    //         shock_results = shock_pro.get_output();

    //         // Clean up after ourselves
    //         shock_pro.clear_input();
    //         shock_pro.clear_output();
            
    //     }
        
  
    // }

    // dbsk2d_shock_storage_sptr shock_storage;
    // shock_storage.vertical_cast(shock_results[0]);
  
    // dbsk2d_ishock_grouping_transform grouper(shock_storage
    //                                          ->get_ishock_graph());
    // grouper.grow_regions();


    // vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_node*> >
    //     fragments = grouper.get_outer_shock_nodes();
    // vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_edge*> >
    //     frag_edges = grouper.get_region_nodes();
    // vcl_map<unsigned int, vcl_vector<dbsk2d_ishock_belm*> > 
    //     frag_belms = grouper.get_region_belms();

    // vgl_polygon<double> model_poly;
    // vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_edge*> >::iterator it;
    // for ( it = frag_edges.begin() ; it != frag_edges.end() ; ++it)
    // {

    //     bool closed_region=(fragments[(*it).first].size()==0)?
    //         true:false;

    //     if ( closed_region && grouper.region_within_image((*it).first,-1))
    //     {
    //         // Grab polygon fragment
    //         grouper.polygon_fragment((*it).first,model_poly);

    //         break;
    //     }
    // }

    // vcl_vector<double> region_stats;
    // vcl_vector<double> app_stats;

    // // get polygon stats
    // grouper.get_region_stats((*it).first,
    //                          model_poly,region_stats);

    // this->get_appearance_stats(frag_edges[(*it).first],
    //                            frag_belms[(*it).first],
    //                            vgl_area(model_poly),
    //                            app_stats);

    // vcl_vector<double> total_stats;
    // total_stats.push_back(0.0); //depth , look at this again
    // total_stats.push_back(1.0); //path prob
    // total_stats.push_back(1.0); //region gap cost
    
    // for ( unsigned int p=0; p < region_stats.size() ; ++p)
    // {
    //     total_stats.push_back(region_stats[p]);
    // }
    
    // for ( unsigned int a=0; a < app_stats.size() ; ++a)
    // {
    //     total_stats.push_back(app_stats[a]);
    // }

    // this->write_output_region_stats(total_stats);

}

// : get closest point
dbsk2d_ishock_bpoint* dbsk2d_transform_manager::get_anchor_pt(
    vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bline*>& pair )
{

    dbsk2d_ishock_bpoint* anchor_pt(0);

    dbsk2d_ishock_bpoint* bp1 = pair.first;
    dbsk2d_ishock_bline*  bl1 = pair.second;

    double d1=vgl_distance(bp1->pt(),bl1->s_pt()->pt());
    double d2=vgl_distance(bp1->pt(),bl1->e_pt()->pt());
    
    double angle1=_vPointPoint(bp1->pt(),bl1->s_pt()->pt());
    double angle2=_vPointPoint(bp1->pt(),bl1->e_pt()->pt());

    vgl_vector_2d<double> vec(vcl_sin(bp1->tangent()),-vcl_cos(bp1->tangent()));
    vgl_line_2d<double> line(bp1->pt(),vec);
    vgl_box_2d<double> box(vgl_point_2d<double>(0,0),vgl_point_2d<double>
                           (image_->ni(),image_->nj()));
    vgl_point_2d<double> lstart,lend;
    bool flag=vgl_intersection(box,line,lstart,lend);

    // convert the pts into bnd_vertex and put into a list
    vcl_vector<vgl_point_2d<double> > bv_list;

    vcl_set<int> local_map;

    vcl_map<int,dbsk2d_ishock_node*> outer_wavefront;
    vcl_map<int,dbsk2d_ishock_node*> source_node;

    // Loop over shock map
    bnd_ishock_map_iter curS = bp1->shock_map().begin();
    for (; curS!=bp1->shock_map().end(); ++curS)
    {
        dbsk2d_ishock_elm* selm = curS->second;
        dbsk2d_ishock_edge* cur_edge = (dbsk2d_ishock_edge*)selm; 
        local_map.insert(cur_edge->rBElement()->id());
        local_map.insert(cur_edge->lBElement()->id());

        if ( cur_edge->cSNode())
        {
            outer_wavefront[cur_edge->cSNode()->id()]=
                cur_edge->cSNode();
        }
        
        if ( cur_edge->pSNode())
        {
            outer_wavefront[cur_edge->pSNode()->id()]=
                cur_edge->pSNode();
        }
        
        if ( cur_edge->pSNode())
        {
            if ( cur_edge->pSNode()->is_a_source())
            {
                source_node[cur_edge->pSNode()->id()]=cur_edge->pSNode();
            }
        }

        if ( cur_edge->cSNode())
        {
            if ( cur_edge->cSNode()->is_a_source())
            {
                source_node[cur_edge->cSNode()->id()]=cur_edge->cSNode();
            }
        }
    }

    

    vcl_map<int,dbsk2d_ishock_node*>::iterator it;
    while ( outer_wavefront.size() > 0 )
    {
        it = outer_wavefront.begin();
        ishock_edge_list pshocks = (*it).second->pShocks();
        dbsk2d_ishock_edge* cshock = (*it).second->cShock();
        dbsk2d_ishock_edge* cshock2 = (*it).second->cShock2();
        
        //remove this edge from the nodes' parent list
        ishock_edge_list::iterator curS = pshocks.begin();
        for(; curS!=pshocks.end(); ++curS)
        {
            dbsk2d_ishock_edge* shock = (*curS);


            if ( local_map.count(shock->lBElement()->id()))
            {
                local_map.insert(shock->rBElement()->id());
            }
           
            if ( local_map.count(shock->rBElement()->id()))
            {
                local_map.insert(shock->lBElement()->id());
            }     
        }

        if ( cshock )
        {
            if ( cshock->cSNode())
            {
                outer_wavefront[cshock->cSNode()->id()]=cshock->cSNode();
            }
            
            dbsk2d_ishock_edge* cur_edge = cshock;
            local_map.insert(cur_edge->rBElement()->id());
            local_map.insert(cur_edge->lBElement()->id());
            
        }


        if ( cshock2 )
        {
            if ( cshock2->cSNode())
            {
                outer_wavefront[cshock2->cSNode()->id()]=cshock2->cSNode();
            }
            
            dbsk2d_ishock_edge* cur_edge = cshock2;
            local_map.insert(cur_edge->rBElement()->id());
            local_map.insert(cur_edge->lBElement()->id());
            
        }

        outer_wavefront.erase(it);
        
    }

    dbsk2d_ishock_bpoint* s_pt=bl1->s_pt();
    dbsk2d_ishock_bpoint* e_pt=bl1->e_pt();

    dbsk2d_ishock_bline* s_pt_line=(
        s_pt->getElmToTheLeftOf(bl1)->id() == bl1->twinLine()->id())?
        (dbsk2d_ishock_bline*)s_pt->getElmToTheRightOf(bl1):
        (dbsk2d_ishock_bline*)s_pt->getElmToTheLeftOf(bl1);

    dbsk2d_ishock_bline* e_pt_line=(
        e_pt->getElmToTheLeftOf(bl1)->id() == bl1->twinLine()->id())?
        (dbsk2d_ishock_bline*)e_pt->getElmToTheRightOf(bl1):
        (dbsk2d_ishock_bline*)e_pt->getElmToTheLeftOf(bl1);


    if ( local_map.count(s_pt_line->id()) && 
         !local_map.count(e_pt_line->id()))
    {
        anchor_pt=bl1->s_pt();
    }
    else if ( local_map.count(e_pt_line->id()) && 
              !local_map.count(s_pt_line->id()))
    {
        anchor_pt=bl1->e_pt();
    }
    else
    {

        bool s_pt_above = _isPointAboveLine(
            bl1->s_pt()->pt(),lstart,lend);
        bool e_pt_above = _isPointAboveLine(
            bl1->e_pt()->pt(),lstart,lend);
        bool source_pt_above = 
            _isPointAboveLine((*source_node.begin()).second->origin(),
                              lstart,lend);
        
        if ( source_pt_above == s_pt_above &&
             source_pt_above != e_pt_above )
        {
            anchor_pt=bl1->s_pt();
        }
        else if ( source_pt_above == e_pt_above &&
                  source_pt_above != s_pt_above )
        {
            anchor_pt=bl1->e_pt();
        }
        else  // This should never happen
        {
            if ( d1 < d2 )
            {
                anchor_pt=bl1->s_pt();
            }
            else
            {
                anchor_pt=bl1->e_pt();
            }
        }
            
    }
 

    return anchor_pt;

}

// chi squared distance
double dbsk2d_transform_manager::min_distance_hist(
    vcl_vector<vgl_point_2d<double> >& foreground,
    vcl_vector<vgl_point_2d<double> >& background,
    vil_image_view<double>& channel,
    double min, double max,unsigned int nbins)
{

    bsta_histogram<double> foreground_pdf(min,max,nbins);
    bsta_histogram<double> background_pdf(min,max,nbins);

    {
        for ( unsigned int f=0; f < foreground.size() ; ++f)
        {
            double x=foreground[f].x();
            double y=foreground[f].y();
            
            double value = vil_bilin_interp_safe_extend(channel,
                                                        x,
                                                        y);
            foreground_pdf.upcount(value,1.0);
        }
    }
   
    {
        for ( unsigned int f=0; f < background.size() ; ++f)
        {
            double x=background[f].x();
            double y=background[f].y();

            double value = vil_bilin_interp_safe_extend(channel,
                                                        x,
                                                        y);
            background_pdf.upcount(value,1.0);
        }
    }

    vcl_vector<double> tempA=foreground_pdf.count_array();
    double norm1=foreground_pdf.area();

    vcl_vector<double> tempB=background_pdf.count_array();
    double norm2=background_pdf.area();
    

    double run_sum=0.0;

    for ( int i=0; i < tempA.size() ; ++i)
    {
        tempA[i]=tempA[i]/norm1;
        tempB[i]=tempB[i]/norm2;
    }

    return min_distance(tempA,tempB);
}

// chi squared distance
double dbsk2d_transform_manager::chi_squared_color_distance(
    vcl_vector<vgl_point_2d<double> >& foreground,
    vcl_vector<vgl_point_2d<double> >& background,
    vil_image_view<double>& channel,
    double min, double max,unsigned int nbins,bool flip)
{

    bsta_histogram<double> foreground_pdf(min,max,nbins);
    bsta_histogram<double> background_pdf(min,max,nbins);

    {
        for ( unsigned int f=0; f < foreground.size() ; ++f)
        {
            double x=foreground[f].x();
            double y=foreground[f].y();
            
            if ( flip )
            {
                x=foreground[f].y();
                y=foreground[f].x();
            }

            double value = vil_bilin_interp_safe_extend(channel,
                                                        x,
                                                        y);
            foreground_pdf.upcount(value,1.0);
        }
    }
   
    {
        for ( unsigned int f=0; f < background.size() ; ++f)
        {
            double x=background[f].x();
            double y=background[f].y();
            
            if ( flip )
            {
                x=background[f].y();
                y=background[f].x();
            }

            double value = vil_bilin_interp_safe_extend(channel,
                                                        x,
                                                        y);
            background_pdf.upcount(value,1.0);
        }
    }

    double chi_squared=0.5*vnl_chi_squared_statistic_12(
        &(foreground_pdf.count_array()[0]),
        &(background_pdf.count_array()[0]),
        nbins,
        true);

    return chi_squared;
}



// chi squared distance

double dbsk2d_transform_manager::mean_LAB_distance(
    vcl_vector<vgl_point_2d<double> >& foreground,
    vcl_vector<vgl_point_2d<double> >& background)
{

    vnl_vector_fixed<double,3> foreground_mean(0.0);
    vnl_vector_fixed<double,3> background_mean(0.0);
    {
        for ( unsigned int f=0; f < foreground.size() ; ++f)
        {
            double x=foreground[f].x();
            double y=foreground[f].y();

            vnl_vector_fixed<double, 3> v;
            v[0] = vil_bilin_interp_safe(L_img_,x,y);
            v[1] = vil_bilin_interp_safe(a_img_,x,y);
            v[2] = vil_bilin_interp_safe(b_img_,x,y);
            foreground_mean += v;
        }
    }

    {
        for ( unsigned int f=0; f < background.size() ; ++f)
        {
            double x=background[f].x();
            double y=background[f].y();

            vnl_vector_fixed<double, 3> v;
            v[0] = vil_bilin_interp_safe(L_img_,x,y);
            v[1] = vil_bilin_interp_safe(a_img_,x,y);
            v[2] = vil_bilin_interp_safe(b_img_,x,y);
            background_mean += v;
        }
    }


    foreground_mean /= foreground.size();
    background_mean /= background.size();

    return distance_LAB(foreground_mean,background_mean,14.0);
}

bool dbsk2d_transform_manager::gap_endpoint(dbsk2d_ishock_bpoint* bp)
{
 
    if (bp->nLinkedElms() <= 2.0)
    {
        return false;
    }

    bool negative=false;
    bool positive=false;
    bool flag=false;
    belm_list::iterator curB = bp->LinkedBElmList.begin();
    for(; curB!=bp->LinkedBElmList.end(); ++curB) 
    {
        dbsk2d_ishock_bline* bline=(dbsk2d_ishock_bline*)(*curB);
        
        if ( bline->get_contour_id() < 0 )
        {
            negative = true;
        }
        else
        {
            positive = true;
        }

    }
    
    if ( positive && negative )
    {
        flag=true;
    }
  
    return flag;


}

// chi squared distance
void dbsk2d_transform_manager::ellipse_fitting(
    vcl_vector<vgl_point_2d<double> >& foreground,
    vcl_vector<double>& ellipse_stats)
{
    double xo=0.0;
    double yo=0.0;
    
    double npts=foreground.size();

    for (unsigned int i=0; i < foreground.size() ; ++i)
    {
        vgl_point_2d<double> pt=foreground[i];

        xo+=pt.x();
        yo+=pt.y();
    }

    double xcenter=xo/npts;
    double ycenter=yo/npts;

    double X2=0.0;
    double Y2=0.0;
    double XY=0.0;

    for ( unsigned int i=0; i < foreground.size() ; ++i)
    {
        vgl_point_2d<double> pt=foreground[i];

        double diffx=pt.x()-xcenter;
        double diffy=pt.y()-ycenter;

        XY += diffx*diffy;
        X2 += diffx*diffx;
        Y2 += diffy*diffy;
    }

    X2=X2/npts;
    Y2=Y2/npts;
    XY=XY/npts;
 
    vnl_double_2x2 Si;
    Si(0,0)=X2;
    Si(0,1)=XY;
    Si(1,0)=XY;
    Si(1,1)=Y2;

    vnl_svd<double> svd(Si);
    
    double major_axis_radius = 2.0*vcl_sqrt(vcl_fabs(svd.W(0)));
    double minor_axis_radius = 2.0*vcl_sqrt(vcl_fabs(svd.W(1)));

    vnl_matrix<double> v=svd.V();

    double orientation(0.0);

    double angle=vcl_fabs(vcl_atan2(v(1,0),v(0,0)));
    
    if ( angle > vnl_math::pi_over_2 )
    {
        orientation=vcl_fmod(vnl_math::pi,angle);
    }
    else
    {
        orientation=angle;
    }

    ellipse_stats.push_back(major_axis_radius*2.0);
    ellipse_stats.push_back(minor_axis_radius*2.0);
    ellipse_stats.push_back(orientation);

    
}


void dbsk2d_transform_manager::grid_points(
    vcl_vector<dbsk2d_ishock_edge*>& region,
    vcl_vector<dbsk2d_ishock_belm*>& belms,
    vcl_vector<vgl_point_2d<double> >& foreground_grid,
    vcl_vector<vgl_point_2d<double> >& background_grid)
{
 
    // Make a dummy coarse shock graph
    dbsk2d_shock_graph_sptr coarse_graph;
    dbsk2d_sample_ishock sampler(coarse_graph);
    sampler.set_sample_resolution(0.5);
    double step_size=1.0;
    
    vcl_set<int> foreground_shocks;

    // Deal with foreground first
    {
        for ( unsigned int i=0; i < region.size() ; ++i)
        {
            dbsk2d_ishock_edge* cur_iedge=region[i];
            
            foreground_shocks.insert(cur_iedge->id());

            // Create a dummy xshock edge
            dbsk2d_shock_node_sptr parent_node = new dbsk2d_shock_node();
            dbsk2d_shock_node_sptr child_node  = new dbsk2d_shock_node();
            dbsk2d_xshock_edge cur_edge(1,parent_node,child_node);
       
            switch (cur_iedge->type())
            {
            case dbsk2d_ishock_elm::POINTPOINT:
                sampler.sample_ishock_edge((dbsk2d_ishock_pointpoint*)
                                           cur_iedge, 
                                           &cur_edge);
                break;
            case dbsk2d_ishock_elm::POINTLINE:
                sampler.sample_ishock_edge((dbsk2d_ishock_pointline*)
                                           cur_iedge, 
                                           &cur_edge);
                break;
            case dbsk2d_ishock_elm::LINELINE:
                sampler.sample_ishock_edge((dbsk2d_ishock_lineline*)
                                           cur_iedge, 
                                           &cur_edge);
                break;
            default:
                break;
            }

            for ( unsigned int s=0; s < cur_edge.num_samples() ; ++s)
            {
                dbsk2d_xshock_sample_sptr sample=cur_edge.sample(s);

                double R1=sample->radius;
                vgl_point_2d<double> pt =sample->pt;
                double theta=sample->theta;
                double phi=0.0;
                double r=step_size;
            

                if (sample->speed != 0 && sample->speed < 99990)
                {
                    phi=vcl_acos(-1.0/sample->speed);
                }
                else
                {
                    phi=vnl_math::pi/2;
                }

                double vec1=theta+phi;
                double vec2=theta-phi;

                while ( r < R1)
                {

                    vgl_point_2d<double> plus_pt=_translatePoint(pt,vec1,r);
                    vgl_point_2d<double> minus_pt=_translatePoint(pt,vec2,r);
                


                    foreground_grid.push_back(plus_pt);
                    foreground_grid.push_back(minus_pt);



                    r+=step_size;
                }
            
                foreground_grid.push_back(pt);

            }
        }

    }

    // Deal with Background second
    {
        // get outer shocks
        vcl_map<int,dbsk2d_ishock_edge*> outer_shocks;

        for ( unsigned int b=0; b < belms.size() ; ++b)
        {
            dbsk2d_ishock_bline* bline=(dbsk2d_ishock_bline*)belms[b];
            dbsk2d_ishock_bline* twinline=bline->twinLine();
            dbsk2d_ishock_bpoint* s_pt=twinline->s_pt();
            dbsk2d_ishock_bpoint* e_pt=twinline->e_pt();
        
            // Twinline first
            {
                bnd_ishock_map shocks=twinline->shock_map();
                bnd_ishock_map_iter curS = shocks.begin();
                for (; curS!=shocks.end(); ++curS)
                {
                    if ( !foreground_shocks.count(curS->second->id()))
                    {
                        outer_shocks[curS->second->id()]=curS->second;
                    } 
                }
            }

            // s_pt second
            {
                bnd_ishock_map shocks=s_pt->shock_map();
                bnd_ishock_map_iter curS = shocks.begin();
                for (; curS!=shocks.end(); ++curS)
                {
                    if ( !foreground_shocks.count(curS->second->id()))
                    {
                        outer_shocks[curS->second->id()]=curS->second;
                    } 
                }
            }

            // e_pt second
            {
                bnd_ishock_map shocks=e_pt->shock_map();
                bnd_ishock_map_iter curS = shocks.begin();
                for (; curS!=shocks.end(); ++curS)
                {
                    if ( !foreground_shocks.count(curS->second->id()))
                    {
                        outer_shocks[curS->second->id()]=curS->second;
                    } 
                }
            }

        }


        vcl_map<int,dbsk2d_ishock_edge*>::iterator it;
        for ( it = outer_shocks.begin() ; it != outer_shocks.end() ; ++it)
        {
            dbsk2d_ishock_edge* cur_iedge=(*it).second;
 
            // Create a dummy xshock edge
            dbsk2d_shock_node_sptr parent_node = new dbsk2d_shock_node();
            dbsk2d_shock_node_sptr child_node  = new dbsk2d_shock_node();
            dbsk2d_xshock_edge cur_edge(1,parent_node,child_node);
       
            switch (cur_iedge->type())
            {
            case dbsk2d_ishock_elm::POINTPOINT:
                sampler.sample_ishock_edge((dbsk2d_ishock_pointpoint*)
                                           cur_iedge, 
                                           &cur_edge);
                break;
            case dbsk2d_ishock_elm::POINTLINE:
                sampler.sample_ishock_edge((dbsk2d_ishock_pointline*)
                                           cur_iedge, 
                                           &cur_edge);
                break;
            case dbsk2d_ishock_elm::LINELINE:
                sampler.sample_ishock_edge((dbsk2d_ishock_lineline*)
                                           cur_iedge, 
                                           &cur_edge);
                break;
            default:
                break;
            }

            for ( unsigned int s=0; s < cur_edge.num_samples() ; ++s)
            {
                dbsk2d_xshock_sample_sptr sample=cur_edge.sample(s);

                double R1=sample->radius;
                vgl_point_2d<double> pt =sample->pt;
                double theta=sample->theta;
                double phi=0.0;
                double r=step_size;
            

                if (sample->speed != 0 && sample->speed < 99990)
                {
                    phi=vcl_acos(-1.0/sample->speed);
                }
                else
                {
                    phi=vnl_math::pi/2;
                }

                double vec1=theta+phi;
                double vec2=theta-phi;

                while ( r < R1)
                {

                    vgl_point_2d<double> plus_pt=_translatePoint(pt,vec1,r);
                    vgl_point_2d<double> minus_pt=_translatePoint(pt,vec2,r);
                



                    if ( plus_pt.x() >= 0 && 
                         plus_pt.y() >= 0 &&
                         plus_pt.x() <= (image_->ni()-1) && 
                         plus_pt.y() <= (image_->nj()-1))
                    {
                        background_grid.push_back(plus_pt);
                        
                    }

                    if ( minus_pt.x() >= 0 && 
                         minus_pt.y() >= 0 &&
                         minus_pt.x() <= (image_->ni()-1) && 
                         minus_pt.y() <= (image_->nj()-1))
                    {
                        background_grid.push_back(minus_pt);
                        
                    }



                    r+=step_size;
                }
            
                if ( pt.x() >= 0 && 
                     pt.y() >= 0 &&
                     pt.x() <= (image_->ni()-1) && 
                     pt.y() <= (image_->nj()-1))
                {
                    background_grid.push_back(pt);
                }
            }
        }
    }


}

double dbsk2d_transform_manager::color_gradient(dbsol_interp_curve_2d_sptr c,
                                                int region_width)
{
    double distance=0.0;
    double delta=0.3;
    double length_threshold=c->length();
    vcl_vector<double> tangents;
    vcl_vector<vsol_point_2d_sptr> curve_pts;
    dbsol_curve_algs::sample(*c, delta, curve_pts, tangents, 
                             length_threshold);  

    for (unsigned j = 0; j < curve_pts.size(); j++) 
    {  
        vgl_point_2d<double> pt = curve_pts[j]->get_p();
        
        //double normal = angle0To2Pi(tangents[j]+vnl_math::pi/2.0f);
        double normal = tangents[j]+vnl_math::pi/2.0f;
        double sum = region_width;

        vgl_point_2d<double> plus_pt(pt.x() + sum*cos(normal), 
                                     pt.y() + sum*sin(normal));
        vgl_point_2d<double> minus_pt(pt.x() - sum*cos(normal), 
                                      pt.y() - sum*sin(normal));

        vnl_vector_fixed<double,3> plus;
        plus[0] = vil_bilin_interp_safe_extend(L_img_,
                                               plus_pt.x(),
                                               plus_pt.y());
        plus[1] = vil_bilin_interp_safe_extend(a_img_,
                                               plus_pt.x(),
                                               plus_pt.y());
        plus[2] = vil_bilin_interp_safe_extend(b_img_,
                                               plus_pt.x(),
                                               plus_pt.y());


        vnl_vector_fixed<double,3> minus;
        minus[0] = vil_bilin_interp_safe_extend(L_img_,
                                                minus_pt.x(),
                                                minus_pt.y());
        minus[1] = vil_bilin_interp_safe_extend(a_img_,
                                                minus_pt.x(),
                                                minus_pt.y());
        minus[2] = vil_bilin_interp_safe_extend(b_img_,
                                                minus_pt.x(),
                                                minus_pt.y());
        
        
       
        vnl_vector_fixed<double, 3> sub = plus-minus;
        double E = sub.two_norm();
        
        distance+=E;
    }

    return distance/c->length();

}

double dbsk2d_transform_manager::likelihood(
    vcl_vector<vgl_point_2d<double> >& curve)
{

    vcl_vector<vsol_point_2d_sptr> pts;

    for ( int i=0; i < curve.size() ; ++i)
    {
        vsol_point_2d_sptr pt=new vsol_point_2d(curve[i]);
        pts.push_back(pt );
    }

    vcl_vector<vsol_point_2d_sptr > region_pts;                 
    dbsol_interp_curve_2d_sptr c = new dbsol_interp_curve_2d();
    dbsol_curve_algs::interpolate_linear(c.ptr(), pts, false);
    dbsol_curve_algs::sample_region_along_curve(*c, 
                                                region_pts, 
                                                0.3f, 
                                                c->length(), 
                                                5.0, 
                                                false);

    vcl_vector<vgl_point_2d<double> > foreground_grid;
    vcl_vector<vgl_point_2d<double> > background_grid;

    for(unsigned i = 0; i < region_pts.size()/2; ++i)
    {
        foreground_grid.push_back(region_pts[i]->get_p());
    }

    for(unsigned i = region_pts.size()/2; i < region_pts.size(); ++i)
    {
        background_grid.push_back(region_pts[i]->get_p());
    }
    

    // 1) Get L difference in a LAB color space
    double L_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               L_img_,
                                               0.0,
                                               100.0,
                                               50);

    // 2) Get a difference in a LAB color space
    double a_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               a_img_,
                                               -110.0,
                                               110.0,
                                               100);
    
    // 3) Get b difference in a LAB color space
    double b_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               b_img_,
                                               -110.0,
                                               110.0,
                                               100);

    // 4) Get b difference in a LAB color space
    double texton_chi2 = chi_squared_color_distance(foreground_grid,
                                                    background_grid,
                                                    texton_image_,
                                                    1.0,
                                                    64.0,
                                                    64);


    double weight1(-0.9521);
    double weight2(-0.6998);
    double weight3(2.7862);
    double weight4(0.6521);
    double weight5(-0.8183);

    double modulus=1*weight1+L_chi2*weight2+a_chi2*weight3+b_chi2*weight4+
        texton_chi2*weight5;

    double sigmoid=1/(1+vcl_exp(-1.0*modulus));

    return sigmoid;

}

double dbsk2d_transform_manager::likelihood(
    vsol_polyline_2d_sptr& curve)
{

    vcl_vector<vsol_point_2d_sptr> pts;

    for ( int i=0; i < curve->size() ; ++i)
    {
        pts.push_back(curve->vertex(i));
    }

    vcl_vector<vsol_point_2d_sptr > region_pts;                 
    dbsol_interp_curve_2d_sptr c = new dbsol_interp_curve_2d();
    dbsol_curve_algs::interpolate_linear(c.ptr(), pts, false);
    dbsol_curve_algs::sample_region_along_curve(*c, 
                                                region_pts, 
                                                0.3f, 
                                                c->length(), 
                                                5.0, 
                                                false);

    vcl_vector<vgl_point_2d<double> > foreground_grid;
    vcl_vector<vgl_point_2d<double> > background_grid;

    for(unsigned i = 0; i < region_pts.size()/2; ++i)
    {
        foreground_grid.push_back(region_pts[i]->get_p());
    }

    for(unsigned i = region_pts.size()/2; i < region_pts.size(); ++i)
    {
        background_grid.push_back(region_pts[i]->get_p());
    }
    

    // 1) Get L difference in a LAB color space
    double L_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               L_img_,
                                               0.0,
                                               1.0,
                                               32);

    // 2) Get a difference in a LAB color space
    double a_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               a_img_,
                                               0.0,
                                               1.0,
                                               32);
    
    // 3) Get b difference in a LAB color space
    double b_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               b_img_,
                                               0.0,
                                               1.0,
                                               32);

    // 4) Get b difference in a LAB color space
    double texton_chi2 = chi_squared_color_distance(foreground_grid,
                                                    background_grid,
                                                    texton_image_,
                                                    1.0,
                                                    64.0,
                                                    64);



    double weight1(-0.9521);
    double weight2(-0.6998);
    double weight3(2.7862);
    double weight4(0.6521);
    double weight5(-0.8183);

    double modulus=1*weight1+L_chi2*weight2+a_chi2*weight3+b_chi2*weight4+
        texton_chi2*weight5;

    double sigmoid=1/(1+vcl_exp(-1.0*modulus));

    return sigmoid;

}

vcl_vector<double> dbsk2d_transform_manager::joint_histogram(
    vcl_vector<vgl_point_2d<double> >& samples,
    vil_image_view<double>& channel1,
    vil_image_view<double>& channel2,
    vil_image_view<double>& channel3,
    double min, double max,unsigned int nbins)
{

    bsta_histogram<double> chan1_hist(min,max,nbins);
    bsta_histogram<double> chan2_hist(min,max,nbins);
    bsta_histogram<double> chan3_hist(min,max,nbins);

    vcl_vector<double> count_array;

    // Per channel
    {
        for ( unsigned int f=0; f < samples.size() ; ++f)
        {
            double y=samples[f].y();
            double x=samples[f].x();
            
            double value1 = vil_bilin_interp_safe_extend(channel1,
                                                        x,
                                                        y);

            double value2 = vil_bilin_interp_safe_extend(channel2,
                                                        x,
                                                        y);

            double value3 = vil_bilin_interp_safe_extend(channel3,
                                                        x,
                                                        y);

            chan1_hist.upcount(value1,1.0);
            chan2_hist.upcount(value2,1.0);
            chan3_hist.upcount(value3,1.0);
        }
    }



    vcl_vector<double> tempA=chan1_hist.count_array();
    vcl_vector<double> tempB=chan2_hist.count_array();
    vcl_vector<double> tempC=chan3_hist.count_array();
   
    double norm=chan1_hist.area()+chan2_hist.area()+chan3_hist.area();
    vcl_vector<double> joint_hist(3*tempA.size(),0.0);

    double run_sum=0.0;

    for ( int i=0; i < tempA.size() ; ++i)
    {
        joint_hist[i]=tempA[i]/norm;
        joint_hist[i+tempA.size()]=tempB[i]/norm;
        joint_hist[i+tempA.size()*2]=tempC[i]/norm;
    }

    return joint_hist;

    

}

double dbsk2d_transform_manager::region_similarity(
    vcl_vector<vgl_point_2d<double> >& 
    region_A_samples,
    vcl_vector<vgl_point_2d<double> >&
    region_B_samples,
    bool color_flag)
{

    // // Debug write out samples
    // {
    //     vcl_ofstream stream("samples_a.txt");
    //     for ( int i=0; i < region_A_samples.size() ; ++i)
    //     {

    //         stream<<region_A_samples[i].x()<<" "<<region_A_samples[i].y()
    //               <<vcl_endl;
    //     }
    //     stream.close();

    // }

    // {
    //     vcl_ofstream stream("samples_b.txt");
    //     for ( int i=0; i < region_B_samples.size() ; ++i)
    //     {

    //         stream<<region_B_samples[i].x()<<" "<<region_B_samples[i].y()
    //               <<vcl_endl;
    //     }
    //     stream.close();

    // }

    double lab_diff(0.0);
    double rgb_diff(0.0);

    if ( color_flag )
    {
    
        // 1) Get concatenated LAB histogram for Region A
        vcl_vector<double> lab_region_A = joint_histogram(region_A_samples,
                                                          L_img_,
                                                          a_img_,
                                                          b_img_,
                                                          0.0,
                                                          1.0,
                                                          25);

        // 2) Get concatenated LAB histogram for Region B
        vcl_vector<double> lab_region_B = joint_histogram(region_B_samples,
                                                          L_img_,
                                                          a_img_,
                                                          b_img_,
                                                          0.0,
                                                          1.0,
                                                          25);

        // 3) Get region diff
        // lab_diff = 0.5*vnl_chi_squared_statistic_12(
        //     &(lab_region_A[0]),
        //     &(lab_region_B[0]),
        //     lab_region_A.size(),
        //     true);
        
        lab_diff = min_distance(lab_region_A,lab_region_B);


        // {
        //     vcl_ofstream stream("hist_a.txt");
        //     for ( int i=0; i < lab_region_A.size() ; ++i)
        //     {
        //         stream<<lab_region_A[i]<<vcl_endl;
        //     }
        //     stream.close();
            
        // }
        
        // {
        //     vcl_ofstream stream("hist_b.txt");
        //     for ( int i=0; i < lab_region_B.size() ; ++i)
        //     {
        //         stream<<lab_region_B[i]<<vcl_endl;
        //     }
        //     stream.close();
            
        // }

    }
    else
    {

        vil_image_resource_sptr img_r = vil_plane(image_, 0);
        vil_image_resource_sptr img_g = vil_plane(image_, 1);
        vil_image_resource_sptr img_b = vil_plane(image_, 2);

        vil_image_view<vxl_byte> red = img_r->get_view();
        vil_image_view<vxl_byte> green = img_g->get_view();
        vil_image_view<vxl_byte> blue = img_b->get_view();

        vil_image_view<double> img_rv;
        vil_image_view<double> img_gv;
        vil_image_view<double> img_bv;

        vil_convert_cast(red,img_rv);
        vil_convert_cast(green,img_gv);
        vil_convert_cast(blue,img_bv);

        // 1) Get concatenated RGB histogram for Region A
        vcl_vector<double> rgb_region_A = joint_histogram(region_A_samples,
                                                          img_rv,
                                                          img_gv,
                                                          img_bv,
                                                          0.0,
                                                          255.0,
                                                          25);

        
        
        // 1) Get concatenated RGB histogram for Region A
        vcl_vector<double> rgb_region_B = joint_histogram(region_B_samples,
                                                          img_rv,
                                                          img_gv,
                                                          img_bv,
                                                          0.0,
                                                          255.0,
                                                          25);
        
        // 3) Get region diff
        rgb_diff = 0.5*vnl_chi_squared_statistic_12(
            &(rgb_region_A[0]),
            &(rgb_region_B[0]),
            rgb_region_A.size(),
            true);

    }


    // 3) Get texton difference in a LAB color space
    double texton_diff = min_distance_hist(region_A_samples,
                                           region_B_samples,
                                           texton_image_,
                                           1.0,
                                           64.0,
                                           64);
    
   
    
    double lambda=0.6;
    double prob(0.0);

    if (color_flag )
    {
        prob=lambda*lab_diff+(1.0-lambda)*texton_diff;
    }
    else
    {
        prob=lambda*rgb_diff+(1.0-lambda)*texton_diff;
    }
    return prob;
}

double dbsk2d_transform_manager::transform_probability(
    vcl_vector<vgl_point_2d<double> >& 
    foreground_grid,
    vcl_vector<vgl_point_2d<double> >&
    background_grid)
{

    // 1) Get L difference in a LAB color space
    double L_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               L_img_,
                                               0.0,
                                               100.0,
                                               50);

    // 2) Get a difference in a LAB color space
    double a_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               a_img_,
                                               -110.0,
                                               110.0,
                                               100);
    
    // 3) Get b difference in a LAB color space
    double b_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               b_img_,
                                               -110.0,
                                               110.0,
                                               100);

    // 4) Get b difference in a LAB color space
    double texton_chi2 = chi_squared_color_distance(foreground_grid,
                                                    background_grid,
                                                    texton_image_,
                                                    1.0,
                                                    64.0,
                                                    64);



    double weight1(-0.9521);
    double weight2(-0.6998);
    double weight3(2.7862);
    double weight4(0.6521);
    double weight5(-0.8183);

    // double weight1(-1.6117);
    // double weight2(-0.0637);
    // double weight3(1.9599);
    // double weight4(2.3581);
    // double weight5(-1.9080);

    // double weight1(-4.5016);
    // double weight2(1.6921);
    // double weight3(0.9562);
    // double weight4(1.0045);
    // double weight5(2.8116);

    double modulus=1*weight1+L_chi2*weight2+a_chi2*weight3+b_chi2*weight4+
        texton_chi2*weight5;

    double sigmoid=1/(1+vcl_exp(-1.0*modulus));

    // vcl_cout<<L_chi2<<" "<<a_chi2<<" "<<b_chi2<<" "<<texton_chi2<<" "<<
    //     sigmoid<<vcl_endl;
    return sigmoid;
}
double dbsk2d_transform_manager::transform_probability(
    double gamma_norm, double k0_norm,double length)
{

    if ( dist_volume_ == 0 )
    {
        return 0.0;
    }

    // Convert
    double gamma_converted  = (gamma_norm+15.0)*2.0;
    double k0_converted     = (k0_norm+15.0)*2.0;
    double length_converted = length*2.0;

    double distance = vil3d_trilin_interp_safe(k0_converted,
                                               gamma_converted,
                                               length_converted,
                                               dist_volume_.origin_ptr(),
                                               dist_volume_.ni(),
                                               dist_volume_.nj(),
                                               dist_volume_.nk(),
                                               dist_volume_.istep(),
                                               dist_volume_.jstep(),
                                               dist_volume_.kstep());
 
    double prob = 1.0-(1.0/(1.0+vcl_exp(distance*logistic_beta0_+
                                    logistic_beta1_)));
    
    return prob;
}

double dbsk2d_transform_manager::transform_probability(
    vcl_vector<vgl_point_2d<double> >& input,
    bool use_length)
{

    // resample curve to always 100 samples

    vcl_vector<vsol_point_2d_sptr> pts;

    for ( int i=0; i < input.size() ; ++i)
    {
        vsol_point_2d_sptr pt=new vsol_point_2d(input[i]);
        pts.push_back(pt );
    }

    dbsol_interp_curve_2d_sptr c = new dbsol_interp_curve_2d();
    dbsol_curve_algs::interpolate_linear(c.ptr(), pts, false);
    dbsol_curve_algs::sample(*c,100,pts);

    if ( gPb_image_ == 0 || c->size() == 0 )
    {
        return 0.0;
    }

    double summation=0.0;
    for ( unsigned int c=0; c < pts.size() ; ++c)
    {
        double x=pts[c]->get_p().x();
        double y=pts[c]->get_p().y();

        double gPb = vil_bilin_interp_safe_extend(gPb_image_,
                                                  x,
                                                  y);

        summation = summation + gPb;
    }
    
    double average_gPb = vcl_min(summation/pts.size(),1.0);

    double sigmoid = average_gPb/normalization_;
    // double length = c->length()/diag_;
    // double sigmoid=0.0;

    // if ( use_length )
    // {
    //     double weight1(-3.0310);
    //     double weight2(10.2010);
    //     double weight3(0.6022);
    //     double modulus=1*weight1+average_gPb*weight2+length*weight3;
    //     sigmoid=1/(1+vcl_exp(-1.0*modulus));

    // }
    // else
    // {
    //     double weight1(-3.0470);
    //     double weight2(10.5959);
    //     double modulus=1*weight1+average_gPb*weight2;
    //     sigmoid=1/(1+vcl_exp(-1.0*modulus));

    // }

    return sigmoid;
}

double dbsk2d_transform_manager::transform_probability(
    vsol_polyline_2d_sptr& input,
    bool use_length)
{
    // resample curve to always 100 samples

    vcl_vector<vsol_point_2d_sptr> pts;

    for ( int i=0; i < input->size() ; ++i)
    {
        pts.push_back(input->vertex(i) );
    }

    dbsol_interp_curve_2d_sptr c = new dbsol_interp_curve_2d();
    dbsol_curve_algs::interpolate_linear(c.ptr(), pts, false);
    dbsol_curve_algs::sample(*c,100,pts);

    if ( gPb_image_ == 0 || c->size() == 0 )
    {
        return 0.0;
    }

    double summation=0.0;
    for ( unsigned int c=0; c < pts.size() ; ++c)
    {
        double x=pts[c]->get_p().x();
        double y=pts[c]->get_p().y();

        double gPb = vil_bilin_interp_safe_extend(gPb_image_,
                                                  x,
                                                  y);

        summation = summation + gPb;
    }
    
    double average_gPb = vcl_min(summation/pts.size(),1.0);
    double length = c->length()/diag_;
    double sigmoid=0.0;

    if ( use_length )
    {
        double weight1(-3.0310);
        double weight2(10.2010);
        double weight3(0.6022);
        double modulus=1*weight1+average_gPb*weight2+length*weight3;
        sigmoid=1/(1+vcl_exp(-1.0*modulus));

    }
    else
    {
        double weight1(-3.0470);
        double weight2(10.5959);
        double modulus=1*weight1+average_gPb*weight2;
        sigmoid=1/(1+vcl_exp(-1.0*modulus));

    }
    return sigmoid;

}

void dbsk2d_transform_manager::start_binary_file(vcl_string binary_file_output)
{
    output_binary_file_ = binary_file_output;

    vcl_ofstream output_binary_file;
    output_binary_file.open(output_binary_file_.c_str(),
                            vcl_ios::out | 
                            vcl_ios::app | 
                            vcl_ios::binary);

    double size_x = image_->ni();
    double size_y = image_->nj();

    output_binary_file.write(reinterpret_cast<char *>(&size_x),
                              sizeof(double));
    output_binary_file.write(reinterpret_cast<char *>(&size_y),
                              sizeof(double));

    output_binary_file.close();


}

void dbsk2d_transform_manager::start_region_file(vcl_string binary_file_output)
{
    output_region_file_ = binary_file_output;

    vcl_ofstream output_region_file;
    output_region_file.open(output_region_file_.c_str(),
                            vcl_ios::out | 
                            vcl_ios::app | 
                            vcl_ios::binary);

    double size_x = image_->ni();
    double size_y = image_->nj();

    output_region_file.write(reinterpret_cast<char *>(&size_x),
                             sizeof(double));
    output_region_file.write(reinterpret_cast<char *>(&size_y),
                             sizeof(double));

    output_region_file.close();


}

void dbsk2d_transform_manager::write_output_polygon(vgl_polygon<double>& poly)
{


    vcl_ofstream output_binary_file;
    output_binary_file.open(output_binary_file_.c_str(),
                            vcl_ios::out | 
                            vcl_ios::app | 
                            vcl_ios::binary);
    
    double num_sheets=poly.num_sheets();
    output_binary_file.write(reinterpret_cast<char *>(&num_sheets),
                              sizeof(double));

    for (unsigned int s = 0; s < poly.num_sheets(); ++s)
    { 

        double num_vertices= poly[s].size();
        output_binary_file.write(reinterpret_cast<char *>(&num_vertices),
                                 sizeof(double));

        for (unsigned int p = 0; p < poly[s].size(); ++p)
        {
            double xcoord=poly[s][p].x();
            double ycoord=poly[s][p].y();

            output_binary_file.write(reinterpret_cast<char *>(&xcoord),
                                     sizeof(double));
            output_binary_file.write(reinterpret_cast<char *>(&ycoord),
                                     sizeof(double));

        }
    }

    output_binary_file.close();


}

void dbsk2d_transform_manager::write_output_region_stats(
    vcl_vector<double>& region_stats)
{


    vcl_ofstream output_binary_file;
    output_binary_file.open(output_region_stats_file_.c_str(),
                            vcl_ios::out | 
                            vcl_ios::app | 
                            vcl_ios::binary);
    

    for ( unsigned int f=0; f < region_stats.size() ; ++f)
    {

        double feature=region_stats[f];
        output_binary_file.write(reinterpret_cast<char *>(&feature),
                                 sizeof(double));

    }
    
    
    output_binary_file.close();


}


void dbsk2d_transform_manager::write_output_region(vgl_polygon<double>& poly)
{


    vcl_ofstream output_region_file;
    output_region_file.open(output_region_file_.c_str(),
                            vcl_ios::out | 
                            vcl_ios::app | 
                            vcl_ios::binary);
    
    double num_contours= (poly[0].size()-1)*4.0+(poly[0].size()-1);
    output_region_file.write(reinterpret_cast<char *>(&num_contours),
                              sizeof(double));
    for (unsigned int p = 0; p < (poly[0].size()-1); ++p)
    {
        double x1_coord = poly[0][p].x();
        double y1_coord = poly[0][p].y();

        double x2_coord = poly[0][p+1].x();
        double y2_coord = poly[0][p+1].y();

        output_region_file.write(reinterpret_cast<char *>(&x1_coord),
                                  sizeof(double));
        output_region_file.write(reinterpret_cast<char *>(&y1_coord),
                                  sizeof(double));

        output_region_file.write(reinterpret_cast<char *>(&x2_coord),
                                  sizeof(double));
        output_region_file.write(reinterpret_cast<char *>(&y2_coord),
                                  sizeof(double));

        double contour_id = 1.0;
        output_region_file.write(reinterpret_cast<char *>(&contour_id),
                                 sizeof(double));
        
    
    }

    output_region_file.close();


}

void dbsk2d_transform_manager::get_appearance_stats(
    vcl_vector<dbsk2d_ishock_edge*>& region,
    vcl_vector<dbsk2d_ishock_belm*>& belms,
    double area,
    vcl_vector<double>& app_stats)
{
    vcl_vector<vgl_point_2d<double> > foreground_grid;
    vcl_vector<vgl_point_2d<double> > background_grid;
    
    this->grid_points(region,
                      belms,
                      foreground_grid,
                      background_grid);
    
    // 1) get region gpb value
    double region_gpb=this->region_gpb_value(foreground_grid)/area;

    // 2) get contour gpb value
    double contour_gpb=this->contour_gpb_value(belms);

    // 3) Get L difference in a LAB color space
    double L_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               L_img_,
                                               0.0,
                                               100.0,
                                               50);

    // 4) Get a difference in a LAB color space
    double a_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               a_img_,
                                               -110.0,
                                               110.0,
                                               100);
    
    // 5) Get b difference in a LAB color space
    double b_chi2 = chi_squared_color_distance(foreground_grid,
                                               background_grid,
                                               b_img_,
                                               -110.0,
                                               110.0,
                                               100);

    // 5) Get b difference in a LAB color space
    double texton_chi2 = chi_squared_color_distance(foreground_grid,
                                                    background_grid,
                                                    texton_image_,
                                                    1.0,
                                                    64.0,
                                                    64);

    // 6) Get mean LAB difference
    double mean_LAB = mean_LAB_distance(foreground_grid,
                                        background_grid);

    // 7) Do ellipse fitting
    vcl_vector<double> ellipse_stats;
    ellipse_fitting(foreground_grid,
                    ellipse_stats);

    app_stats.push_back(region_gpb);
    app_stats.push_back(contour_gpb);
    app_stats.push_back(L_chi2);
    app_stats.push_back(a_chi2);
    app_stats.push_back(b_chi2);
    app_stats.push_back(texton_chi2);
    app_stats.push_back(mean_LAB);

    for ( unsigned int i=0; i < ellipse_stats.size() ; ++i)
    {
        app_stats.push_back(ellipse_stats[i]); // major axis length
    }

}

void dbsk2d_transform_manager::get_extra_belms(
    vcl_vector<dbsk2d_ishock_belm*>& region_belms,
    vcl_set<int>& key,
    vcl_set<int>& closed_region_key,
    vcl_map<int,dbsk2d_ishock_bline*>& output_lines)
{
    vcl_vector<dbsk2d_ishock_belm*>::iterator lit;  
    for (lit = region_belms.begin() ; lit != region_belms.end() ; ++lit)
    {
        dbsk2d_ishock_bline* bline = (dbsk2d_ishock_bline*)(*lit);
        dbsk2d_ishock_bpoint* s_pt = bline->s_pt();
        dbsk2d_ishock_bpoint* e_pt = bline->e_pt();
        output_lines[bline->id()]=bline;
        closed_region_key.insert(bline->id());

        if ( bline->get_contour_id() >= 0 )
        {
            key.insert(bline->id());
        }

        if ( gap_endpoint(s_pt) )
        {
            key.insert(s_pt->id());
        }

        if ( gap_endpoint(e_pt) )
        {
            key.insert(e_pt->id());
        }

    }     

    for (lit = region_belms.begin() ; lit != region_belms.end() ; ++lit)
    {
    
        dbsk2d_ishock_bline* bline = (dbsk2d_ishock_bline*)(*lit);

        {
            dbsk2d_ishock_bpoint* s_pt=bline->s_pt();

            dbsk2d_ishock_bline* left_element_spt =
                (dbsk2d_ishock_bline*)(s_pt->getElmToTheRightOf(bline));
            dbsk2d_ishock_bline* right_element_spt =
                (dbsk2d_ishock_bline*)(s_pt->getElmToTheLeftOf(bline));

            if ( !output_lines.count(left_element_spt->twinLine()->id()))
            { 
                output_lines[left_element_spt->id()]=left_element_spt;
            }
           
            if ( !output_lines.count(right_element_spt->twinLine()->id()))
            { 
                output_lines[right_element_spt->id()]=right_element_spt;
            }
            
        }

        
        {
            dbsk2d_ishock_bpoint* e_pt=bline->e_pt();
           
            dbsk2d_ishock_bline* left_element_ept =
                (dbsk2d_ishock_bline*)(e_pt->getElmToTheRightOf(bline));
            dbsk2d_ishock_bline* right_element_ept =
                (dbsk2d_ishock_bline*)(e_pt->getElmToTheLeftOf(bline));

            if ( !output_lines.count(left_element_ept->twinLine()->id()))
            { 
                output_lines[left_element_ept->id()]=left_element_ept;
            }
           
            if ( !output_lines.count(right_element_ept->twinLine()->id()))
            { 
                output_lines[right_element_ept->id()]=right_element_ept;
            }
        }
   
    }

    // vcl_map<int, dbsk2d_ishock_bline*>::iterator oit;
    // for (oit = output_lines.begin() ; oit != output_lines.end() ; ++oit)
    // {
    //     key.insert((*oit).first);
    // }
}

void dbsk2d_transform_manager::write_output_region(
    vcl_vector<dbsk2d_ishock_belm*>& region_belms)
{


    vcl_ofstream output_region_file;
    output_region_file.open(output_region_file_.c_str(),
                            vcl_ios::out | 
                            vcl_ios::app | 
                            vcl_ios::binary);
    
    vcl_map<int, dbsk2d_ishock_bline*> output_lines;
  
    vcl_vector<dbsk2d_ishock_belm*>::iterator lit;  
    for (lit = region_belms.begin() ; lit != region_belms.end() ; ++lit)
    {
        dbsk2d_ishock_bline* bline = (dbsk2d_ishock_bline*)(*lit);
        output_lines[bline->id()]=bline;
    }     

    // for (lit = region_belms.begin() ; lit != region_belms.end() ; ++lit)
    // {
    
    //     dbsk2d_ishock_bline* bline = (dbsk2d_ishock_bline*)(*lit);

    //     {
    //         dbsk2d_ishock_bpoint* s_pt=bline->s_pt();
    //         belm_list linkedbelms=s_pt->LinkedBElmList;
            
    //         belm_list::iterator bit;
    //         for ( bit = linkedbelms.begin() ; bit != linkedbelms.end() ; ++bit)
    //         {
    //             dbsk2d_ishock_bline* bline_extra = (dbsk2d_ishock_bline*)(*bit);
    //             if ( !output_lines.count(bline_extra->twinLine()->id()))
    //             { 
    //                 output_lines[bline_extra->id()]=bline_extra;
    //             }
                
    //         }
    //     }

        
    //     {
    //         dbsk2d_ishock_bpoint* e_pt=bline->e_pt();
    //         belm_list linkedbelms=e_pt->LinkedBElmList;
            
    //         belm_list::iterator bit;
    //         for ( bit = linkedbelms.begin() ; bit != linkedbelms.end() ; ++bit)
    //         {
    //             dbsk2d_ishock_bline* bline_extra = (dbsk2d_ishock_bline*)(*bit);
    //             if ( !output_lines.count(bline_extra->twinLine()->id()))
    //             { 
    //                 output_lines[bline_extra->id()]=bline_extra;
    //             }
    //         }
    //     }
   
    // }

        
 
    double num_contours= output_lines.size()*4.0+output_lines.size();
    output_region_file.write(reinterpret_cast<char *>(&num_contours),
                              sizeof(double));
  
    vcl_map<int, dbsk2d_ishock_bline*>::iterator oit;
    for (oit = output_lines.begin() ; oit != output_lines.end() ; ++oit)
    {

        dbsk2d_ishock_bline* bline = (dbsk2d_ishock_bline*)((*oit).second);
                
        dbsk2d_ishock_bpoint* source = bline->s_pt();
        dbsk2d_ishock_bpoint* target = bline->e_pt();

        double x1_coord = source->pt().x();
        double y1_coord = source->pt().y();

        double x2_coord = target->pt().x();
        double y2_coord = target->pt().y();

        output_region_file.write(reinterpret_cast<char *>(&x1_coord),
                                  sizeof(double));
        output_region_file.write(reinterpret_cast<char *>(&y1_coord),
                                  sizeof(double));

        output_region_file.write(reinterpret_cast<char *>(&x2_coord),
                                  sizeof(double));
        output_region_file.write(reinterpret_cast<char *>(&y2_coord),
                                  sizeof(double));

        double contour_id = (*oit).second->get_contour_id();
        output_region_file.write(reinterpret_cast<char *>(&contour_id),
                                 sizeof(double));

    }
    output_region_file.close();

}


void dbsk2d_transform_manager::write_output_region(
    vcl_vector<dbsk2d_bnd_contour_sptr>& contours,
    vcl_vector<vgl_point_2d<double> >& gap_filler)
{


    vgl_fit_lines_2d<double> fitter;
    fitter.set_min_fit_length(2);
    fitter.set_rms_error_tol(0.05f);
    fitter.add_curve(gap_filler);
    fitter.fit();

    vcl_vector<vgl_line_segment_2d<double> > segs;
    segs= fitter.get_line_segs();
    
    vcl_ofstream output_region_file;
    output_region_file.open(output_region_file_.c_str(),
                            vcl_ios::out | 
                            vcl_ios::app | 
                            vcl_ios::binary);
    
    vcl_map<int,dbsk2d_bnd_edge_sptr> output_edges;
    double contour_id_orig=contours[0]->get_id();

    for ( unsigned int c=0; c < contours.size() ; ++c)
    {
        dbsk2d_bnd_contour_sptr curve = contours[c];
        for ( unsigned int v=0; v < curve->num_edges() ; ++v)
        {
            dbsk2d_bnd_edge_sptr edge=curve->bnd_edge(v);
            output_edges[edge->get_id()]=edge;
        }
    }
        
 
    double num_contours= output_edges.size()*4.0+output_edges.size()+
        segs.size()*4.0 +segs.size();
    output_region_file.write(reinterpret_cast<char *>(&num_contours),
                              sizeof(double));
  
    vcl_map<int, dbsk2d_bnd_edge_sptr>::iterator oit;
    for (oit = output_edges.begin() ; oit != output_edges.end() ; ++oit)
    {

        dbsk2d_bnd_edge_sptr bedge = (*oit).second;
                
        dbsk2d_ishock_bpoint* source = bedge->bnd_v1()->bpoint();
        dbsk2d_ishock_bpoint* target = bedge->bnd_v2()->bpoint();

        double x1_coord = source->pt().x();
        double y1_coord = source->pt().y();

        double x2_coord = target->pt().x();
        double y2_coord = target->pt().y();

        output_region_file.write(reinterpret_cast<char *>(&x1_coord),
                                  sizeof(double));
        output_region_file.write(reinterpret_cast<char *>(&y1_coord),
                                  sizeof(double));

        output_region_file.write(reinterpret_cast<char *>(&x2_coord),
                                  sizeof(double));
        output_region_file.write(reinterpret_cast<char *>(&y2_coord),
                                  sizeof(double));

        double contour_id = (*oit).second->get_id();
        output_region_file.write(reinterpret_cast<char *>(&contour_id),
                                 sizeof(double));

    }

    for ( unsigned int ls=0; ls < segs.size() ; ++ls)
    {

        vgl_line_segment_2d<double> segment=segs[ls];

        double x1_coord = segment.point1().x();
        double y1_coord = segment.point1().y();

        double x2_coord = segment.point2().x();
        double y2_coord = segment.point2().y();

        output_region_file.write(reinterpret_cast<char *>(&x1_coord),
                                  sizeof(double));
        output_region_file.write(reinterpret_cast<char *>(&y1_coord),
                                  sizeof(double));

        output_region_file.write(reinterpret_cast<char *>(&x2_coord),
                                  sizeof(double));
        output_region_file.write(reinterpret_cast<char *>(&y2_coord),
                                  sizeof(double));

        output_region_file.write(reinterpret_cast<char *>(&contour_id_orig),
                                 sizeof(double));

    }
    output_region_file.close();

}

void dbsk2d_transform_manager::write_output_polygon(
    vcl_vector<dbsk2d_bnd_contour_sptr>& contours,
    vcl_vector<vgl_point_2d<double> >& gap_filler)
{


   
    vcl_vector<vgl_point_2d<double> > hull_points;

    for (unsigned int k=0; k < gap_filler.size() ; ++k)
    {

        hull_points.push_back(gap_filler[k]);
    }

    for ( unsigned int c=0; c < contours.size() ; ++c)
    {
        dbsk2d_bnd_contour_sptr curve = contours[c];
        for ( unsigned int v=0; v < curve->num_edges()+1 ; ++v)
        {
            dbsk2d_bnd_vertex_sptr vertex=curve->bnd_vertex(v);
            hull_points.push_back(vertex->point());
        }
    }

    vgl_convex_hull_2d<double> convex_hull(hull_points);
    vgl_polygon<double> poly=convex_hull.hull();

    vcl_ofstream output_binary_file;
    output_binary_file.open(output_binary_file_.c_str(),
                            vcl_ios::out | 
                            vcl_ios::app | 
                            vcl_ios::binary);
    
    double num_vertices= poly[0].size();
    output_binary_file.write(reinterpret_cast<char *>(&num_vertices),
                              sizeof(double));

    for (unsigned int p = 0; p < poly[0].size(); ++p)
    {
        double xcoord = poly[0][p].x();
        double ycoord = poly[0][p].y();

        output_binary_file.write(reinterpret_cast<char *>(&xcoord),
                                  sizeof(double));
        output_binary_file.write(reinterpret_cast<char *>(&ycoord),
                                  sizeof(double));
        
    
    }

    output_binary_file.close();

}


void dbsk2d_transform_manager::save_image_poly(
    vgl_polygon<double>& vgl_poly,
    vcl_string filename)
{
    
    vil_image_resource_sptr img_r = vil_plane(image_, 0);
    vil_image_resource_sptr img_g = vil_plane(image_, 1);
    vil_image_resource_sptr img_b = vil_plane(image_, 2);

    vsol_polygon_2d_sptr vsol_poly = bsol_algs::poly_from_vgl(vgl_poly);
    
    vsol_poly->compute_bounding_box();
    vsol_box_2d_sptr bbox = vsol_poly->get_bounding_box();
    double minx = bbox->get_min_x()-5 < 0 ? 0 : bbox->get_min_x()-5;
    double miny = bbox->get_min_y()-5 < 0 ? 0 : bbox->get_min_y()-5;

    vil_image_view<vil_rgb<vxl_byte> > 
        temp((int)vcl_ceil(bbox->width() + 10), 
             (int)vcl_ceil(bbox->height() + 10), 1); 
    vil_rgb<vxl_byte> bg_col(255, 255, 255);
    temp.fill(bg_col);

    vil_image_view<vxl_byte> img_rv = img_r->get_view();
    vil_image_view<vxl_byte> img_gv = img_g->get_view();
    vil_image_view<vxl_byte> img_bv = img_b->get_view();

    // do not include boundary
    vgl_polygon_scan_iterator<double> psi(vgl_poly, false);  
    for (psi.reset(); psi.next(); ) 
    {
        int y = psi.scany();
        for (int x = psi.startx(); x <= psi.endx(); ++x) 
        {
            if (x < 0 || y < 0)
            {
                continue;
            }
            if (x >= int(img_r->ni()) || y >= int(img_r->nj()))
            { 
                continue;
            }
            int xx = (int)vcl_floor(x - minx + 0.5); 
            int yy = (int)vcl_floor(y - miny + 0.5);
            if (xx < 0 || yy < 0)
            {
                continue;
            }
            if (double(xx) > bbox->width() || double(yy) > bbox->height())
            {
                continue;
            }
            temp(xx,yy) = 
                vil_rgb<vxl_byte>(img_rv(x,y), img_gv(x,y), img_bv(x,y));
        }
    }

    vil_image_resource_sptr out_img = vil_new_image_resource_of_view(temp);
    vil_save_image_resource(out_img, 
                            filename.c_str()); 

}
