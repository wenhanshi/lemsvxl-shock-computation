#include "dbnl_eno.h"
#include <vcl_cstdio.h>
#include <vcl_cassert.h>
#include <vcl_limits.h>
#include <vcl_algorithm.h>  

//:
// \file
// \author Ricardo Fabbri

#if !VCL_STATIC_CONST_INIT_FLOAT_NO_DEFN
const double  dbnl_eno::near_zero_value
      VCL_STATIC_CONST_INIT_FLOAT_DEFN(1.0e-8);
const double  dbnl_eno::interval_tolerance
      VCL_STATIC_CONST_INIT_FLOAT_DEFN(1.0e-8);
#endif

//:
// \return true on success
bool dbnl_eno_1d::
interpolate(double const *data, unsigned len, ptrdiff_t stride)
{
   /* if don't have enough data to interpolate intervals, nothing to
      do here */
#ifndef NDEBUG
   if (len < DBNL_ENO_DATA_LENGTH) {
      // better err treatment TODO
      printf("data gridline too short to interpolate.\n");
      printf("in: enointerv_compute_interpolants\n");
      return false;
   }
#endif


   double nhood[DBNL_ENO_DATA_LENGTH]; /* storage for copying
                                           local neighborhood */

   vcl_vector<dbnl_eno_interp>::iterator interp_ptr;

   stride_ = stride;
   interp_.resize(len*stride -1);
   interp_ptr = interp_.begin();

   /* intialize loop index before calculations */
   unsigned ii = 0;

   /* A macro for improved readability */
#define ENOStoreDataList4(l,f0,f1,f2,f3) \
    do{(l)[0] = (f0); (l)[1] = (f1); (l)[2] = (f2); (l)[3] = (f3);}while(0)


   /* do the first border value by hand */
   ENOStoreDataList4(nhood,
                     border_value_, data[0],
                     data[stride * 1], data[stride * 2]);
   interp_ptr->set_bounds(ii, ii + 1);
   interp_ptr->interpolate(nhood);
   data += stride;
   interp_ptr += stride;

   for (ii = 1; ii +2 < len; ++ii) {

      ENOStoreDataList4(nhood,
                        data[stride * -1], data[0],
                        data[stride * 1], data[stride * 2]);
      interp_ptr->set_bounds(ii, ii + 1);
      interp_ptr->interpolate(nhood);

      /* increment to next element */
      data += stride;
      interp_ptr += stride;
   }

   /* do the last border value by hand */
   ENOStoreDataList4(nhood,
                     data[stride * -1], data[0],
                     data[stride * 1], border_value_);
   interp_ptr->set_bounds(ii, ii + 1);
   interp_ptr->interpolate(nhood);

#undef ENOStoreDataList4

   abscissas_.set_size(len);
   for (ii=0; ii< len; ++ii)
      abscissas_[ii] = ii;

   return true;
}

//:
// stride does not apply to abcissa vector x
// \return true on success
bool dbnl_eno_1d::
interpolate(double const *data, double const *x, unsigned len, ptrdiff_t stride)
{
   /* if don't have enough data to interpolate intervals, nothing to
      do here */
#ifndef NDEBUG
   if (len < DBNL_ENO_DATA_LENGTH) {
      // better err treatment TODO
      printf("data gridline too short to interpolate.\n");
      printf("in: enointerv_compute_interpolants\n");
      return false;
   }
#endif

   double nhood[DBNL_ENO_DATA_LENGTH]; /* storage for copying
                                           local neighborhood */

   vcl_vector<dbnl_eno_interp>::iterator interp_ptr;

   stride_ = stride;
   interp_.resize(len*stride -1);
   interp_ptr = interp_.begin();

   /* intialize loop index before calculations */
   unsigned ii = 0;

   /* A macro for improved readability */
#define ENOStoreDataList4(l,f0,f1,f2,f3) \
    do{(l)[0] = (f0); (l)[1] = (f1); (l)[2] = (f2); (l)[3] = (f3);}while(0)


   interp_ptr->set_bounds(x[0], x[1]);
   interp_ptr->interpolate(data[0], data[stride*1], data[stride*2],x[0],x[1],x[2]);

   data += stride;
   interp_ptr += stride;

   for (ii = 1; ii +2 < len; ++ii) {

      ENOStoreDataList4(nhood,
                        data[stride * -1], data[0],
                        data[stride * 1], data[stride * 2]);

      interp_ptr->set_bounds(x[ii], x[ii + 1]);
      interp_ptr->interpolate(nhood,x+ii-1);

      /* increment to next element */
      data += stride;
      interp_ptr += stride;
   }

   interp_ptr->set_bounds(x[ii], x[ii + 1]);
   interp_ptr->interpolate(*data,data[stride*1],data[-1*stride],x[ii],x[ii+1],x[ii-1]);

#undef ENOStoreDataList4

   abscissas_.set_size(len);
   abscissas_.copy_in(x);

   return true;
}

void dbnl_eno_1d::
print(vcl_ostream& strm) const
{

   strm << "==== Eno 1D ====" << vcl_endl
        << "len: " << size() << vcl_endl;

   unsigned i;
   for (i=0; i < size(); ++i)
      interp_.at(i).print(strm);
}

double dbnl_eno_1d::
sample(double x) const
{
   unsigned i=interval_index(x);
   return interp_[i].sample(x);
}

bool dbnl_eno_1d::
sample(unsigned size, vnl_vector<double> &f, vnl_vector<double> &xx) const
{
   f.set_size(size);
   xx.set_size(size);

   double pace= ( x(this->size()) - x(0) ) / (size-1);

   for (unsigned i=0; i<size; i++) {
      xx[i] = x(0) + i*pace;
      f[i] = sample(xx(i));
   }

   return true;
}

//: Finds index of the interval containing a given abscissa, using
// binary search.
// If s correspond to one of the endpoints/samples, we return the index of the
// interval to the "left" of it.
unsigned dbnl_eno_1d::
interval_index(double x) const
{
   // binary search for s in vector of arclens
   const vnl_vector<double>::const_iterator
      p = vcl_lower_bound(abscissas_.begin(), abscissas_.end(), x);

   unsigned i = p - abscissas_.begin();

   return (i == 0)? i : i-1;
}


// copied from dbnl_eno_poly.cxx & dbnl_eno_intrep.cxx
// there is not any instruction about implementation of sample() and interpolate in .h
// here are these implementations

//: interpolation code for non-integral coords
void dbnl_eno_interp ::
interpolate(double const data[DBNL_ENO_DATA_LENGTH], double const x[DBNL_ENO_DATA_LENGTH])
{
   double
           a2,      /* 2nd order coefficient of forward polynomial */
           b2,      /* 2nd order coefficient of backward polynomial */
           d21,d32,d31,
           det;

   /* compute leading coefficient of forward and backward polynomials */
   d21 = x[1]-x[2];
   d32 = x[2]-x[3];
   d31 = x[1]-x[3];

   det= d21*d31*d32;
   a2 = (d32*data[1] -d31*data[2] +d21*data[3]  )/det;

   d21 = x[0]-x[1];
   d32 = x[1]-x[2];
   d31 = x[0]-x[2];

   det= d21*d31*d32;
   b2 = (d32*data[0] -d31*data[1] +d21*data[2]  )/det;

   /* determine which direction to use for interpolation */
   forward_ = vcl_fabs(a2) < vcl_fabs(b2);

   /* choose polynomial with smaller variation, where variation is
      measured as absolute value of leading polynomial coefficient.*/
   a2 = forward_ ? a2:b2;

   /* compute and store all polynomial coefficients for this interpolant */
   coeffs_[dbnl_eno_poly::second_order_index] = a2;

   a2 = -a2;
   d21 = x[2]-x[1];
   coeffs_[dbnl_eno_poly::first_order_index] =
           (data[2] - data[1])/d21 + a2*(x[1] + x[2]);

   coeffs_[dbnl_eno_poly::zero_order_index] =
           data[1] + x[1]*(a2*x[1] - coeffs_[dbnl_eno_poly::first_order_index]);

   order_=2;
}

//: fit poly to 3 pts
void dbnl_eno_interp::
interpolate(double d1, double d2, double d3, double x1, double x2, double x3)
{
   double
           a2,      /* 2nd order coefficient of forward polynomial */
           d21,d32,d31,
           det;

   forward_ = x3 > x1;

   // compute leading coefficient
   d21 = x1-x2;
   d32 = x2-x3;
   d31 = x1-x3;

   det= d21*d31*d32;
   a2 = ( d32*d1 -d31*d2 +d21*d3 )/det;

   // compute and store all polynomial coefficients for this interpolant
   coeffs_[dbnl_eno_poly::second_order_index] = a2;

   a2 = -a2;
   d21 = -d21;
   coeffs_[dbnl_eno_poly::first_order_index] =
           (d2 - d1)/d21 + a2*(x1 + x2);

   coeffs_[dbnl_eno_poly::zero_order_index] =
           d1 + x1*(a2*x1 - coeffs_[dbnl_eno_poly::first_order_index]);

   order_=2;
}

//: specific interpolation code for integer coords
void dbnl_eno_interp ::
interpolate(double const data[DBNL_ENO_DATA_LENGTH])
{
   double
           a2,      /* 2nd order coefficient of forward polynomial */
           b2,      /* 2nd order coefficient of backward polynomial */
           c2;      /* 2nd order coefficient in choosen direction */
   int const
           off = 1;      /* offset in data array */

   /* compute leading coefficient of forward and backward polynomials */
   a2 = (data[off+2] - 2.0* data[off+1] + data[off+0])/2.0;
   b2 = (data[off+1] - 2.0* data[off+0] + data[off-1])/2.0;

   /* determine which direction to use for interpolation */
   forward_ = vcl_fabs(a2) < vcl_fabs(b2);

   /* choose polynomial with smaller variation, where variation is
      measured as absolute value of leading polynomial coefficient.*/
   c2 = forward_ ? a2:b2;

   /* compute and store all polynomial coefficients for this interpolant */
   coeffs_[dbnl_eno_poly::second_order_index] = c2;
   coeffs_[dbnl_eno_poly::first_order_index] = data[off+1] - (c2*(2*start_+1) + data[off+0]);
   coeffs_[dbnl_eno_poly::zero_order_index] =
           (data[off+0] + c2*start_*(start_+1)) - (start_*(data[off+1] - data[off+0]));

   order_=2;
}

void dbnl_eno_interp::
print(vcl_ostream &strm) const
{
   strm << "==== Interpolant ====" << vcl_endl
        << "forward: " << (forward_?"true":"false") << vcl_endl
        << "start,end:  " << start_ << " , " << end_ << vcl_endl;

   dbnl_eno_poly::print(strm);
}

// Constants
#if !VCL_STATIC_CONST_INIT_INT_NO_DEFN
const unsigned dbnl_eno_poly::zero_order_index
VCL_STATIC_CONST_INIT_INT_DEFN(0);
const unsigned  dbnl_eno_poly::first_order_index
VCL_STATIC_CONST_INIT_INT_DEFN(1);
const unsigned  dbnl_eno_poly::second_order_index
VCL_STATIC_CONST_INIT_INT_DEFN(2);
#endif

void dbnl_eno_poly::
print(vcl_ostream &strm) const
{
   strm << "==== Polynomial ====" << vcl_endl
        << "order: " << order() << vcl_endl
        << "coefficients: ";

   for (unsigned i=0; i<=order(); ++i)
      strm << coeff(i) << " ";
   strm << vcl_endl;
}

dbnl_eno_poly operator-(const dbnl_eno_poly &f1, const dbnl_eno_poly &f2)
{
   dbnl_eno_poly diff_poly(2);

   diff_poly[f1.second_order_index] = f1[f1.second_order_index] - f2[f1.second_order_index];
   diff_poly[f1.first_order_index]  = f1[f1.first_order_index] - f2[f1.first_order_index];
   diff_poly[f1.zero_order_index]   = f1[f1.zero_order_index] - f2[f1.zero_order_index];

   return diff_poly;
}

double dbnl_eno_poly::
sample(double x) const
{
   return x*(coeffs_[second_order_index]*x + coeffs_[first_order_index])
          + coeffs_[zero_order_index];
}

// TODO
// - an eval method would be useful for testing...
// - use a class from vnl?
