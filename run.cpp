// SPDX-License-Identifier: MIT
// author: Emmanuel Le Trong <emmanuel.le-trong@cnrs-orleans.fr>
#include <cassert>
#include <vector>
#include <numbers>
    using std::numbers::pi;
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include <fmt/ostream.h>
    using fmt::println;


/* Solves f(x) = 0 for x, which should be bracketed by [x_min, x_max]. Uses 
   the GNU Scientific Library (GSL).
  
   Parameters:
   ==========
    - f        : the function whose root we're after. Its signature should be
                 double f (double x, void* params)
    - params   : parameters to pass to the function.
    - x_min    : lower bound for the root.
    - x_max    : upper bound for the root.
    - max_iter : maximum number of iterations the function will attempt to find
                 a root.
    - rtol     : relative tolerance, see the GSL documentation for gsl_roots.h.
    - atol     : absolute tolerance, see the GSL documentation for gsl_roots.h.
  
   Returns: a tuple with
   =======
    - root   : The root if found, the closest estimate if not. NaN if there is
               not a single root inside the [x_min, x_max] interval.
    - iter   : The number of iterations the function performed.
    - status : An error code indicating if a solution has been found, with
               possible values:
                - GSL_SUCCESS if a root has been found.
                - GSL_EDOM if there is not a single root inside the 
                  [x_min, x_max] interval.
                - Any other value if the procedure otherwise fails (see the GSL 
                  documentation for gsl_root.h).
 */
    auto
fsolve1 (
      auto f
    , auto& params
    , double x_min
    , double x_max
    , int max_iter = 100
    , double rtol = 1e-8
    , double atol = 0.
){
    {
            const auto
        f1 = f (x_min, &params);
            const auto
        f2 = f (x_max, &params);
        if (f1 * f2 > 0) 
        {
                return std::tuple { 
                  std::nan ("")
                , 0
                , static_cast <int> (GSL_EDOM)
            };
        }
    }
        auto
    solver = gsl_root_fsolver_alloc (gsl_root_fsolver_bisection);
        auto
    func = gsl_function {
          .function = f
        , .params = &params
    };
    gsl_root_fsolver_set (solver, &func, x_min, x_max);

        int
    status;
        auto
    iter = 0;
        double
      root
    , x_lo
    , x_hi
    ;
    do 
    {
        ++iter;
        status = gsl_root_fsolver_iterate (solver);
        root = gsl_root_fsolver_root (solver);
        x_lo = gsl_root_fsolver_x_lower (solver);
        x_hi = gsl_root_fsolver_x_upper (solver);
        status = gsl_root_test_interval (x_lo, x_hi, atol, rtol);
    } 
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (solver);
    return std::tuple { root, iter, status };
}

// Tests for the function fsolve1
    namespace
tests
{
        bool
    approximatively_equal (
          double x
        , double y
        , double rtol = 1e-6
        , double atol = 0.
    ){
        return fabs (x - y) < fmin (abs (x), abs (y)) * rtol + atol;
    }
        struct
    params1_t
    {
            double
          a
        , b
        , c
        ;
    };
        double
    f1 (double x, void *params)
    {
            const auto
        p = static_cast <params1_t const*> (params);
        return (p->a * x + p->b) * x + p->c;
    }

        auto
    test_fsolve1 ()
    {
            auto
        params = params1_t { 1., 0., -5. };
            const auto
        [ root, iter, status ] = fsolve1 (f1, params, 0., 5.);
        assert (status == GSL_SUCCESS);
            const auto
        expected = sqrt (5.);
        /*
        println (
             "Found {} in {} iterations, expected {}"
            , root
            , iter
            , expected
        );
        */
        assert (approximatively_equal (root, expected));
    }
} // namespace tests

/* Solves a system of two (possibly non-linear) equations, starting from an
   initial guess (x0, x1). Uses the GNU Scientific Library (GSL).
  
   Parameters:
   ==========
    - f        : the function whose roots we're after. Its signature should be
                    int f (const gsl_vector *x, void *params, gsl_vector *f)
                 where 
                  - x is the current solution estimate (input), 
                  - params are the parameters of f (input),
                  - f are the calculated values of the function at the current 
                    estimate (output)
                 the function shall return GSL_SUCCESS or any other value to 
                 indicate a failure (see the GSL documentation for 
                 gsl_multiroots.h).
    - params   : parameters to pass to the function.
    - x1       : initial guess for the first root
    - x2       : initial guess for the second root
    - max_iter : maximum number of iterations the function will attempt to find
                 a root.
    - atol     : absolute tolerance, see the GSL documentation for 
                 gsl_multiroots.h.
  
   Returns: a tuple with
   =======
    - x1     : The first root if found, the closest estimate if not.
    - x2     : The second root if found, the closest estimate if not.
    - iter   : The number of iterations the function performed.
    - status : An error code indicating if a solution has been found, with
               possible values:
                - GSL_SUCCESS if a root has been found.
                - Any other value if the procedure otherwise fails (see the GSL 
                  documentation for gsl_multiroot.h).
 */
    auto
fsolve2 (
      auto f
    , auto& params
    , double x1
    , double x2
    , int max_iter = 100
    , double atol = 1e-8
){
        auto
    solver = gsl_multiroot_fsolver_alloc (gsl_multiroot_fsolver_hybrids, 2);
        auto
    func = gsl_multiroot_function {
          .f = f
        , .n = 2
        , .params = &params
    };
        auto
    x_init = gsl_vector_alloc (2);
    gsl_vector_set (x_init, 0, x1);
    gsl_vector_set (x_init, 1, x2);
    gsl_multiroot_fsolver_set (solver, &func, x_init);

        int
    status;
        auto
    iter = 0;
    do
    {
        ++iter;
        status = gsl_multiroot_fsolver_iterate (solver);
        if (status) break;
        status = gsl_multiroot_test_residual(solver->f, atol);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
        const auto
    r = std::tuple { 
          gsl_vector_get (solver->x, 0)
        , gsl_vector_get (solver->x, 1)
        , iter
        , status
    };
    gsl_multiroot_fsolver_free (solver);
    gsl_vector_free (x_init);
    return r;
}

// Tests for fsolve2
    namespace
tests
{
        struct
    params2_t
    {
            double
          a
        , b
        ;
    };

        int
    f2 (const gsl_vector *x, void *params, gsl_vector *f)
    {
            const auto
        p = static_cast <params2_t const*> (params);
            const auto
        x0 = gsl_vector_get (x, 0);
            const auto
        x1 = gsl_vector_get (x, 1);
            const auto
        f0 = p->a * (1. - x0);
            const auto
        f1 = p->b * (x1 - x0 * x0);
        gsl_vector_set (f, 0, f0);
        gsl_vector_set (f, 1, f1);
        return GSL_SUCCESS;
    }

        auto
    test_fsolve2 ()
    {
            auto
        params = params2_t { 1., 10. };
            const auto
        [ x0, x1, iter, status] = fsolve2 (f2, params, -10., -5.);
        assert (status == GSL_SUCCESS);
            const auto
        expected0 = 1.;
            const auto
        expected1 = 1.;
        /*
        println (
             "Found ({},{}) in {} iterations, expected ({},{})"
            , x0
            , x1
            , iter
            , expected0
            , expected1
        );
        */
        assert (approximatively_equal (x0, expected0));
        assert (approximatively_equal (x1, expected1));
    }
} // namespace tests

/* Computes the radius on the interface given the radii of the bubble and drop,
   and the surface tensions.

   Parameters:
   ==========
    - R_b : the radius of the bubble,
    - R_m : the radius of the drop,
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension.

   Returns: a tuple with
    - into_bubble : true if the interface intrudes into the bubble, false 
                    otherwise,
    - R_i         : the radius of the interface.
   =======
 */
    auto
calc_R_i (
      double R_b
    , double R_m
    , double sigma_bl
    , double sigma_ml
    , double sigma_bm)
{
        const auto
    r = R_b * sigma_bl - R_m * sigma_ml;
    return std::tuple { r > 0, R_m * R_b * sigma_bm / fabs (r) };
}

/* Computes the volume of a spherical cap.

   Parameters:
   ==========
    - r : the radius of the cap,
    - h : the height of the cap,
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension.

   Returns: the volume of the cap.
   =======
 */
    auto
calc_V_c (double r, double h)
{ 
    return pi * pow (h, 2) * (r - h / 3);
}
/* Computes the surface area of a spherical cap.

   Parameters:
   ==========
    - r : the radius of the cap,
    - h : the height of the cap.
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension.

   Returns: the surface area of the cap.
   =======
 */
    auto
calc_A_c (double r, double h)
{
    return 2 * pi * r * h;
}


// Parameters for the following function.
    struct
f_R_i_params_t 
{
        double
      r_b
    , r_m
    , sigma_bl
    , sigma_ml
    , sigma_bm
    , theta_b
    , theta_m
    ;
};

// The function whose roots are the radii of the bubble and the drop after 
// attachment.
    int
f_R_i (const gsl_vector *x, void *params, gsl_vector *f)
{
        const auto
    p = static_cast <f_R_i_params_t const*> (params);
        const auto
    R_b = gsl_vector_get (x, 0);
        const auto
    R_m = gsl_vector_get (x, 1);
        const auto
    [ is_right, R_i ] = calc_R_i (R_b, R_m, p->sigma_bl, p->sigma_ml, p->sigma_bm);
        const auto
    alpha = p->theta_b + p->theta_m;
        const auto
    a2 = -0.5 * pow (R_b, 2) * pow (R_m, 2) * pow (sin (alpha), 2) 
        / (pow (R_b, 2) + pow (R_m, 2) + 2 * R_b * R_m * cos (alpha))
    ;
        const auto
    h_i_1 = R_i - sqrt (pow (R_i, 2)- a2);
        const auto
    V_c_i = calc_V_c (R_i, h_i_1);

        const auto
    h_b_2 = R_b + sqrt (pow (R_b, 2) - a2);
        const auto
    V_b = calc_V_c (R_b, h_b_2) + (is_right ? -V_c_i : V_c_i);
        const auto
    V0_b = 4. / 3 * pi * pow (p->r_b, 3);

        const auto
    h_m_2 = R_m + sqrt (pow (R_m, 2) - a2);
        const auto
    V_m = calc_V_c (R_m, h_m_2) + (is_right ? V_c_i : -V_c_i);
        const auto
    V0_m = 4. / 3 * pi * pow (p->r_m, 3);

    gsl_vector_set (f, 0, V0_b - V_b);
    gsl_vector_set (f, 1, V0_m - V_m);

    return GSL_SUCCESS;
}

/* Computes the radii of the bubble and drop after attachment, given the radii
   before (and the surface tensions).

   Parameters:
   ==========
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension,
    - R0_b     : radius of the bubble before attachment,
    - R0_m     : radius of the drop before attachment.

   Returns: a tuple with:
   =======
    - theta_b : the contact angle inside the bubble,
    - theta_m : the contact angle inside the drop,
    - R_b     : the radius of the bubble,
    - R_m     : the radius of the drop,
    - iter    : the iteration count required to compute the solution,
    - status  : the error code, see the documentation for fsolve2.
 */
    auto
calc_R_bm (
      double sigma_bl
    , double sigma_ml
    , double sigma_bm
    , double R0_b
    , double R0_m
){
        const auto
    theta_b = acos (
          (pow (sigma_ml, 2) - pow (sigma_bl, 2) - pow (sigma_bm, 2)) 
        / 2 / sigma_bl / sigma_bm
    );
        const auto
    theta_m = acos (
          (pow (sigma_bl, 2) - pow (sigma_bm, 2) - pow (sigma_ml, 2)) 
        / 2 / sigma_bm / sigma_ml
    );
        auto
    params = f_R_i_params_t 
    {
          .r_b = R0_b
        , .r_m = R0_m
        , .sigma_bl = sigma_bl 
        , .sigma_ml = sigma_ml 
        , .sigma_bm = sigma_bm 
        , .theta_b  = theta_b  
        , .theta_m  = theta_m  
    };
        const auto
    [ R_b, R_m, iter, status] = fsolve2 (f_R_i, params, R0_b, R0_m);
    return std::tuple {
          theta_b
        , theta_m
        , R_b
        , R_m
        , iter
        , status
    };
}

/* Computes the forces exerted on a system of bubble and drop.

   Parameters:
   ==========
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension,
    - R0_b     : radius of the bubble before attachment,
    - R0_m     : radius of the drop before attachment,
    - rho_b    : density of the fluid inside the bubble,
    - rho_m    : density of the molten metal inside the drop,
    - rho_l    : density of the surrounding liquid,
    - g        : gravitational acceleration.

   Returns: a tuple with:
   =======
    - F_attach   : the force of attachment,
    - F_buoyancy : the net buoyancy of the system,
    - F_detach   : the force of detachment.
 */
    auto
calc_forces (
      double sigma_bl
    , double sigma_ml
    , double sigma_bm
    , double R0_b
    , double R0_m
    , double rho_b
    , double rho_m
    , double rho_l
    , double g
){
        const auto
    [ theta_b, theta_m, R_b, R_m, iter, status] = calc_R_bm (
          sigma_bl
        , sigma_ml
        , sigma_bm
        , R0_b
        , R0_m
    );
        const auto
    [ is_right, R_i ] = calc_R_i (R_b, R_m, sigma_bl, sigma_ml, sigma_bm);
        const auto
    alpha = theta_b + theta_m;
        const auto
    c = sqrt (pow (R_b, 2) + pow (R_m, 2) + 2 * R_b * R_m * cos (alpha));
        const auto
    a2 = -0.5 * pow (R_b, 2) * pow (R_m, 2) * pow (sin (alpha), 2) 
        / (pow (R_b, 2) + pow (R_m, 2) + 2 * R_b * R_m * cos (alpha))
    ;
        const auto
    h_i_1 = R_i - sqrt (pow (R_i, 2)- a2);
        const auto
    h_b_2 = R_b + sqrt (pow (R_b, 2) - a2);
        const auto
    V_c_i = calc_V_c (R_i, h_i_1);
        const auto
    V_b = calc_V_c (R_b, h_b_2) + (is_right ? -V_c_i : V_c_i);

        const auto
    h_m_2 = R_m + sqrt (pow (R_m, 2) - a2);
        const auto
    V_m = calc_V_c (R_m, h_m_2) + (is_right ? V_c_i : -V_c_i);

        const auto
    A0_b = 4 * pi * pow (R0_b, 2);
        const auto
    A0_m = 4 * pi * pow (R0_m, 2);

        const auto
    A_b = calc_A_c (R_b, h_b_2);
        const auto
    A_m = calc_A_c (R_m, h_m_2);
        const auto
    A_i = calc_A_c (R_i, h_i_1);

        const auto
    DG = (A_b - A0_b) * sigma_bl + (A_m - A0_m) * sigma_ml + A_i * sigma_bm;
        const auto
    Dl = R0_b + R0_m - c;
        const auto
    F_attach = DG / Dl;

        const auto
    F_b = V_b * (rho_l - rho_b) * g;
        const auto
    F_m = V_m * (rho_l - rho_m) * g;
    assert (F_b > 0);
    assert (F_m < 0);

        const auto
    F_buoyancy = F_b + F_m;
        const auto
    F_detach = F_b - F_m;

    return std::tuple { F_attach, F_buoyancy, F_detach };
}

// Parameters for the following function.
    struct
f_r_b_params_t
{
        double
      sigma_lb
    , sigma_lm
    , sigma_bm
    , R0_b
    , R0_m
    , rho_b
    , rho_m
    , rho_l
    , g
    ;
};

// Function whose root is the radius of the bubble before attachment that
// provides zero buoyancy for the bubble-drop system.
    double
f_zero_buoyancy_b (double x, void *params)
{
        const auto
    p = static_cast <f_r_b_params_t const*> (params);
        const auto
    R0_b = x;
        const auto
    [ F_attach, F_buoyancy, F_detach ] = calc_forces (
          p -> sigma_lb
        , p -> sigma_lm
        , p -> sigma_bm
        , R0_b
        , p -> R0_m
        , p -> rho_b
        , p -> rho_m
        , p -> rho_l
        , p -> g
    );
    return F_buoyancy;
}
/* Computes the radius of the bubble before attachment that provides zero 
   buoyancy for the bubble-drop system.

   Parameters:
   ==========
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension,
    - R0_b_lo  : radius of the bubble before attachment, lower bound,
    - R0_b_hi  : radius of the bubble before attachment, higher bound,
    - R0_m     : radius of the drop before attachment,
    - rho_b    : density of the fluid inside the bubble,
    - rho_m    : density of the molten metal inside the drop,
    - rho_l    : density of the surrounding liquid,
    - g        : gravitational acceleration.

   Returns: the tuple returned by fsolve1, see above.
   =======
 */
    auto
calc_R_b_zero_buoyancy (
      double sigma_lb
    , double sigma_lm
    , double sigma_bm
    , double R0_b_lo
    , double R0_b_hi
    , double R0_m
    , double rho_b
    , double rho_m
    , double rho_l
    , double g
){
        auto
    params = f_r_b_params_t
    {
          .sigma_lb = sigma_lb
        , .sigma_lm = sigma_lm
        , .sigma_bm = sigma_bm
        , .R0_b     = 0. // not used
        , .R0_m     = R0_m     
        , .rho_b    = rho_b   
        , .rho_m    = rho_m   
        , .rho_l    = rho_l   
        , .g        = g       
    };
    return fsolve1 (f_zero_buoyancy_b, params, R0_b_lo, R0_b_hi);
}

// Function whose root is the radius of the drop before attachment that
// provides zero buoyancy for the bubble-drop system.
    double
f_zero_buoyancy_m (double x, void *params)
{
        const auto
    p = static_cast <f_r_b_params_t const*> (params);
        const auto
    R0_m = x;
        const auto
    [ F_attach, F_buoyancy, F_detach ] = calc_forces (
          p -> sigma_lb
        , p -> sigma_lm
        , p -> sigma_bm
        , p -> R0_b
        , R0_m
        , p -> rho_b
        , p -> rho_m
        , p -> rho_l
        , p -> g
    );
    return F_buoyancy;
}
/* Computes the radius of the drop before attachment that provides zero 
   buoyancy for the bubble-drop system.

   Parameters:
   ==========
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension,
    - R0_b     : radius of the bubble before attachment,
    - R0_m_lo  : radius of the drop before attachment, lower bound,
    - R0_m_hi  : radius of the drop before attachment, higher bound,
    - rho_b    : density of the fluid inside the bubble,
    - rho_m    : density of the molten metal inside the drop,
    - rho_l    : density of the surrounding liquid,
    - g        : gravitational acceleration.

   Returns: the tuple returned by fsolve1, see above.
   =======
 */
    auto
calc_R_m_zero_buoyancy (
      double sigma_lb
    , double sigma_lm
    , double sigma_bm
    , double R0_b
    , double R0_m_lo
    , double R0_m_hi
    , double rho_b
    , double rho_m
    , double rho_l
    , double g
){
        auto
    params = f_r_b_params_t
    {
          .sigma_lb = sigma_lb
        , .sigma_lm = sigma_lm
        , .sigma_bm = sigma_bm
        , .R0_b     = R0_b
        , .R0_m     = 0. // not used
        , .rho_b    = rho_b   
        , .rho_m    = rho_m   
        , .rho_l    = rho_l   
        , .g        = g       
    };
    return fsolve1 (f_zero_buoyancy_m, params, R0_m_lo, R0_m_hi);
}

// Function whose roots are the radii of the bubble and drop before attachment 
// that provides zero buoyancy for the bubble-drop system at the detachment
// limit.
    int
f_zero_buoyancy_bm (const gsl_vector *x, void *params, gsl_vector *f)
{
        const auto
    p = static_cast <f_r_b_params_t const*> (params);
        const auto
    r_b = gsl_vector_get (x, 0);
        const auto
    r_m = gsl_vector_get (x, 1);
        const auto
    [ F_attach, F_buoyancy, F_detach ] = calc_forces (
          p -> sigma_lb
        , p -> sigma_lm
        , p -> sigma_bm
        , r_b
        , r_m
        , p -> rho_b
        , p -> rho_m
        , p -> rho_l
        , p -> g
    );
    gsl_vector_set (f, 0, F_detach - F_attach);
    gsl_vector_set (f, 1, F_buoyancy);
    return GSL_SUCCESS;
}
/* Computes radii of the bubble and drop before attachment that provides zero 
   buoyancy for the bubble-drop system at the detachment limit.
  
   Parameters:
   ==========
    - sigma_bl  : bubble-liquid surface tension,
    - sigma_ml  : drop-liquid surface tension,
    - sigma_bm  : bubble-drop surface tension,
    - R0_b_init : initial guess for the radius of the bubble before attachment,
    - R0_m_init : initial guess for the radius of the drop before attachment,
    - rho_b     : density of the fluid inside the bubble,
    - rho_m     : density of the molten metal inside the drop,
    - rho_l     : density of the surrounding liquid,
    - g         : gravitational acceleration.

   Returns: the tuple returned by fsolve2, see above.
   =======
 */
    auto
calc_R_bm_zero_buoyancy (
      double sigma_lb
    , double sigma_lm
    , double sigma_bm
    , double R0_b_init
    , double R0_m_init
    , double rho_b
    , double rho_m
    , double rho_l
    , double g
){
        auto
    params = f_r_b_params_t
    {
          .sigma_lb = sigma_lb
        , .sigma_lm = sigma_lm
        , .sigma_bm = sigma_bm
        , .R0_b     = 0. // not used
        , .R0_m     = 0. // not used
        , .rho_b    = rho_b   
        , .rho_m    = rho_m   
        , .rho_l    = rho_l   
        , .g        = g       
    };
    return fsolve2 (f_zero_buoyancy_bm, params, R0_b_init, R0_m_init);
}

// Function whose root is the radius of the bubble before attachment 
// such that the force of detachment balances exactly the force of attachment.
    double
f_att_eq_det_b (double x, void *params)
{
        const auto
    p = static_cast <f_r_b_params_t const*> (params);
        const auto
    r_b = x;
        const auto
    [ F_attach, F_buoyancy, F_detach ] = calc_forces (
          p -> sigma_lb
        , p -> sigma_lm
        , p -> sigma_bm
        , r_b
        , p -> R0_m
        , p -> rho_b
        , p -> rho_m
        , p -> rho_l
        , p -> g
    );
    return F_attach - F_detach;
}
/* Computes radius of the bubble before attachment at the detachment limit.
  
   Parameters:
   ==========
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension,
    - R0_b_lo  : radius of the bubble, lower bound
    - R0_b_hi  : radius of the bubble, upper bound
    - R0_m     : radius of the drop,
    - rho_b    : density of the fluid inside the bubble,
    - rho_m    : density of the molten metal inside the drop,
    - rho_l    : density of the surrounding liquid,
    - g        : gravitational acceleration.

   Returns: the tuple returned by fsolve1, see above.
   =======
 */
    auto
calc_R_b_detach (
      double sigma_lb
    , double sigma_lm
    , double sigma_bm
    , double R0_b_lo
    , double R0_b_hi
    , double R0_m
    , double rho_b
    , double rho_m
    , double rho_l
    , double g
){
        auto
    params = f_r_b_params_t
    {
          .sigma_lb = sigma_lb
        , .sigma_lm = sigma_lm
        , .sigma_bm = sigma_bm
        , .R0_b      = 0. // not used
        , .R0_m      = R0_m     
        , .rho_b    = rho_b   
        , .rho_m    = rho_m   
        , .rho_l    = rho_l   
        , .g        = g       
    };
    return fsolve1 (f_att_eq_det_b, params, R0_b_lo, R0_b_hi);
}
// Function whose root is the radius of the drop before attachment 
// such that the force of detachment balances exactly the force of attachment.
    double
f_att_eq_det_m (double x, void *params)
{
        const auto
    p = static_cast <f_r_b_params_t const*> (params);
        const auto
    R0_m = x;
        const auto
    [ F_attach, F_buoyancy, F_detach ] = calc_forces (
          p -> sigma_lb
        , p -> sigma_lm
        , p -> sigma_bm
        , p -> R0_b
        , R0_m
        , p -> rho_b
        , p -> rho_m
        , p -> rho_l
        , p -> g
    );
    return F_attach - F_detach;
}
/* Computes radius of the bubble before attachment at the detachment limit.
  
   Parameters:
   ==========
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension,
    - R0_b     : radius of the bubble,
    - R0_m_lo  : radius of the drop, lower bound
    - R0_m_hi  : radius of the drop, upper bound
    - rho_b    : density of the fluid inside the bubble,
    - rho_m    : density of the molten metal inside the drop,
    - rho_l    : density of the surrounding liquid,
    - g        : gravitational acceleration.

   Returns: the tuple returned by fsolve1, see above.
   =======
 */
    auto
calc_R_m_detach (
      double sigma_lb
    , double sigma_lm
    , double sigma_bm
    , double R0_b
    , double R0_m_lo
    , double R0_m_hi
    , double rho_b
    , double rho_m
    , double rho_l
    , double g
){
        auto
    params = f_r_b_params_t
    {
          .sigma_lb = sigma_lb
        , .sigma_lm = sigma_lm
        , .sigma_bm = sigma_bm
        , .R0_b     = R0_b
        , .R0_m     = 0. // not used
        , .rho_b    = rho_b   
        , .rho_m    = rho_m   
        , .rho_l    = rho_l   
        , .g        = g       
    };
    return fsolve1 (f_att_eq_det_m, params, R0_m_lo, R0_m_hi);
}

/* Produces the data needed to plot the fig. 6 of the paper.
   Parameters:
   ==========
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension,
    - rho_b    : density of the fluid inside the bubble,
    - rho_m    : density of the molten metal inside the drop,
    - rho_l    : density of the surrounding liquid,
    - g        : gravitational acceleration.
    - R0_b_min : bubble radius range, lower bound
    - R0_b_max : bubble radius range, upper bound
    - count_b  : bubble radius count
    - R0_m_min : drop radius range, lower bound 
    - R0_m_max : drop radius range, upper bound 
    - count_m  : drop radius count

   Returns: nothing.
   =======

   Side effects: generates the following CSV files:
   ============
    - "<name>_zero_buoyancy_b.csv"
    - "<name>_detachment1_b.csv"
    - "<name>_detachment2_m.csv"
   which are used to plot the fig. 6. (See the fig6.gnuplot file.)
 
 */
    void
fig6 (
      std::string const& name
    , double sigma_lb
    , double sigma_lm
    , double sigma_bm
    , double rho_b
    , double rho_m
    , double rho_l
    , double g
    , double R0_b_min
    , double R0_b_max
    , size_t count_b
    , double R0_m_min
    , double R0_m_max
    , size_t count_m
){
        const auto
    delta_b = pow (R0_b_max / R0_b_min, 1. / static_cast <double> (count_b - 1));
        const auto
    delta_m = pow (R0_m_max / R0_m_min, 1. / static_cast <double> (count_m - 1));
    {
            auto
        o_buo = std::ofstream { fmt::format ("fig_6_{}_zero_buoyancy_b.csv", name) };
            auto
        o_det1 = std::ofstream { fmt::format ("fig_6_{}_detachment1_b.csv", name) };
            auto
        tip_is_done = false;
            double
          r_b_tip
        , r_m_tip
        ;
        for (auto i_m = 0u; i_m < count_m; ++i_m)
        {
                const auto
            R0_m = R0_m_min * pow (delta_m, i_m);
                const auto
            [ R0_b_buo, iter_buo, status_buo ] = calc_R_b_zero_buoyancy (
                  sigma_lb
                , sigma_lm
                , sigma_bm
                , R0_b_min
                , R0_b_max
                , R0_m
                , rho_b
                , rho_m
                , rho_l
                , g
            );
            if (status_buo != GSL_SUCCESS)
            {
                // no zero buoyancy in the [r_b_min, r_b_max] range
                continue;
            }
                const auto
            [ F_attach, F_buoyancy, F_detach ] = calc_forces (
                  sigma_lb
                , sigma_lm
                , sigma_bm
                , R0_b_buo
                , R0_m
                , rho_b
                , rho_m
                , rho_l
                , g
            );
            if (F_detach <= F_attach)
            {
                println (o_buo, "{} {}", R0_m, R0_b_buo);
            }
            else if (not tip_is_done)
            {
                    int
                  iter
                , status
                ;
                    std::tie
                (r_b_tip, r_m_tip, iter, status) = calc_R_bm_zero_buoyancy (
                      sigma_lb
                    , sigma_lm
                    , sigma_bm
                    , R0_b_buo // r_b_0
                    , R0_m // r_m_0
                    , rho_b
                    , rho_m
                    , rho_l
                    , g
                );
                assert (status == GSL_SUCCESS);
                println (o_buo, "{} {}", r_m_tip, r_b_tip);
                tip_is_done = true;
            }
                const auto
            [ R0_b_det1, iter_det1, status_det1 ] = calc_R_b_detach (
                  sigma_lb
                , sigma_lm
                , sigma_bm
                , R0_b_buo
                , R0_b_max
                , R0_m
                , rho_b
                , rho_m
                , rho_l
                , g
            );
            if (status_det1 == GSL_SUCCESS) 
            {
                println (o_det1, "{} {}", R0_m, R0_b_det1);
            }
        }
        // Connect to the tip
        println (o_det1, "{} {}", r_m_tip, r_b_tip);
    }{
            auto
        o_det2 = std::ofstream { fmt::format ("fig_6_{}_detachment2_m.csv", name) };
            auto
        tip_is_done = false;
            double
          R0_b_tip
        , R0_m_tip
        ;
        for (auto i_b = 0u; i_b < count_b; ++i_b)
        {
                const auto
            R0_b = R0_b_min * pow (delta_b, i_b);
                const auto
            [ R0_m_buo, iter_buo, status_buo ] = calc_R_m_zero_buoyancy (
                  sigma_lb
                , sigma_lm
                , sigma_bm
                , R0_b
                , R0_m_min
                , R0_m_max
                , rho_b
                , rho_m
                , rho_l
                , g
            );
            if (status_buo != GSL_SUCCESS)
            {
                // no zero buoyancy in the [r_b_min, r_b_max] range
                continue;
            }
                const auto
            [ F_attach, F_buoyancy, F_detach ] = calc_forces (
                  sigma_lb
                , sigma_lm
                , sigma_bm
                , R0_b
                , R0_m_buo 
                , rho_b
                , rho_m
                , rho_l
                , g
            );
            if (F_detach <= F_attach)
            {
                // noop
            }
            else if (not tip_is_done)
            {
                    int
                  iter
                , status
                ;
                    std::tie
                (R0_b_tip, R0_m_tip, iter, status) = calc_R_bm_zero_buoyancy (
                      sigma_lb
                    , sigma_lm
                    , sigma_bm
                    , R0_b // r_b_0
                    , R0_m_buo // r_m_0
                    , rho_b
                    , rho_m
                    , rho_l
                    , g
                );
                assert (status == GSL_SUCCESS);
                tip_is_done = true;
            }
                const auto
            [ R0_m_det2, iter_det2, status_det2 ] = calc_R_m_detach (
                  sigma_lb
                , sigma_lm
                , sigma_bm
                , R0_b
                , R0_m_buo
                , R0_m_max
                , rho_b
                , rho_m
                , rho_l
                , g
            );
            if (status_det2 == GSL_SUCCESS) 
            {
                println (o_det2, "{} {}", R0_m_det2, R0_b);
            }
        }
        // Refine around the tip
        {
                const auto
            [ R0_m_det2, iter_det2, status_det2 ] = calc_R_m_detach (
                  sigma_lb
                , sigma_lm
                , sigma_bm
                , R0_b_tip * 0.999 // Slightly below to actually get a result
                , R0_m_tip
                , R0_m_max
                , rho_b
                , rho_m
                , rho_l
                , g
            );
            assert (status_det2 == GSL_SUCCESS);
            println (o_det2, "{} {}", R0_m_det2, R0_b_tip);
            
        }
        // Connect to the tip
        println (o_det2, "{} {}", R0_m_tip, R0_b_tip);
    }
}

/* Produces the data needed to plot the fig. 3a of the paper.
   Parameters:
   ==========
    - sigma_bl : bubble-liquid surface tension,
    - sigma_ml : drop-liquid surface tension,
    - sigma_bm : bubble-drop surface tension,
    - rho_b    : density of the fluid inside the bubble,
    - rho_m    : density of the molten metal inside the drop,
    - rho_l    : density of the surrounding liquid,
    - g        : gravitational acceleration.
    - R0_b_min : bubble radius range, lower bound
    - R0_b_max : bubble radius range, upper bound
    - count_b  : bubble radius count

   Returns: nothing.
   =======

   Side effects: generates the following CSV files:
   ============
   which are used to plot the fig. 3a. (See the fig3a.gnuplot file.)
 
 */
    void
fig3a (
      std::string const& name
    , double sigma_lb
    , double sigma_lm
    , double sigma_bm
    , double rho_b
    , double rho_m
    , double rho_l
    , double g
    , double R0_b_min
    , double R0_b_max
    , size_t count_b
){
        const auto
    delta_b = pow (R0_b_max / R0_b_min, 1. / static_cast <double> (count_b - 1));
        auto
    o = std::ofstream { fmt::format ("fig_3a_{}.csv", name) };
    for (const auto R0_m: { 4e-6, 16e-6, 64e-6, 256e-6, 1000e-6, 5000e-6 })
    {
        println (o, "# R0_m = {}", R0_m);
        for (auto i_b = 0u; i_b < count_b; ++i_b)
        {
                const auto
            R0_b = R0_b_min * pow (delta_b, i_b);
                const auto
            [ F_attach, F_buoyancy, F_detach ] = calc_forces (
                  sigma_lb
                , sigma_lm
                , sigma_bm
                , R0_b
                , R0_m
                , rho_b
                , rho_m
                , rho_l
                , g
            );
            println (
                  o
                , "{} {}"
                , R0_b
                , F_attach / F_detach
            );
        }
        println (o, "\n");
    }
}
    int
main ()
{
    tests::test_fsolve1 ();
    tests::test_fsolve2 ();

        const auto
    planets = std::array { "Asteroid", "Vesta", "Moon", "Mars", "Earth" };
        const auto
    gs = std::array { 0.02, 0.25, 1., 3., 10. };
        const auto
    rho_bs = std::array { 
          14.976759863428649
        , 14.977963101198736
        , 14.981886697004017
        , 14.99234957710964
        , 15.028969175914598
    };

    for (auto i = 0u; i < size (planets); ++i)
    {
            const auto&
        planet = planets.at (i);
            const auto
        g = gs.at (i);
            const auto
        rho_b = rho_bs.at (i);
        fig6 (
              planet
            , 0.37  // sigma_lb
            , 1.50  // sigma_lm
            , 1.20  // sigma_bm
            , rho_b
            , 7874  // rho_m
            , 2800  // rho_l
            , g
            , 1e-9  // r_b_min
            , 1e0   // r_b_max
            , 1000  // count_b
            , 1e-7  // r_m_min
            , 1e-1  // r_m_max
            , 700   // count_m
        );
        fig3a (
              planet
            , 0.37  // sigma_lb
            , 1.50  // sigma_lm
            , 1.20  // sigma_bm
            , rho_b
            , 7874  // rho_m
            , 2800  // rho_l
            , g
            , 0.4e-6  // r_b_min
            , 1e-2   // r_b_max
            , 100  // count_b
        );
    }
    return 0;
}
