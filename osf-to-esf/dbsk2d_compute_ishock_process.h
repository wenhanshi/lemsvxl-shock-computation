//
// Created by wenhan on 7/19/17.
//

#ifndef DBSK2D_ISHOCK_COMPUTE_DBSK2D_ISHOCK_COMPUTER_H
#define DBSK2D_ISHOCK_COMPUTE_DBSK2D_ISHOCK_COMPUTER_H


#include "bpro1/bpro1_process.h"
#include "bpro1/bpro1_parameters.h"
#include <vsol/vsol_polygon_2d_sptr.h>
#include <vsol/vsol_polyline_2d_sptr.h>

class dbsk2d_compute_ishock_process : public bpro1_process
{
public:

    dbsk2d_compute_ishock_process();
    virtual ~dbsk2d_compute_ishock_process();
    //: Clone the process
    virtual bpro1_process* clone() const;

    vcl_string name();

    vcl_vector< vcl_string > get_input_type();
    vcl_vector< vcl_string > get_output_type();

    int input_frames();
    int output_frames();

    bool execute();
    bool finish();

    vsol_polygon_2d_sptr smooth_closed_contour(vsol_polygon_2d_sptr polygon);

    vsol_polyline_2d_sptr add_noise_to_contour(vsol_polyline_2d_sptr poly, double noise_radius);
    vsol_polygon_2d_sptr add_noise_to_contour(vsol_polygon_2d_sptr poly, double noise_radius);
    vsol_polyline_2d_sptr fit_lines_to_contour(vsol_polyline_2d_sptr poly, double rms);
    vsol_polygon_2d_sptr fit_lines_to_contour(vsol_polygon_2d_sptr poly, double rms);

protected:

private:
};


#endif //DBSK2D_ISHOCK_COMPUTE_DBSK2D_ISHOCK_COMPUTER_H
