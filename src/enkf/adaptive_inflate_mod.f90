module adaptive_inflate_mod
! module: adaptive_inflate    Adaptive inflation 
!                             modified from DART
!
! prgmmr: jjh                 date: 2017-02-20  
!> Operations and storage required for various adaptive inflation algorithms

use constants, only: pi
use kinds, only: r_kind, i_kind, r_double, r_single
use params, only: prior_inf_damping, prior_inf_sd_lower_bound, prior_inf_lower_bound, prior_inf_upper_bound

implicit none
private

public :: update_inflation,          solve_quadratic

contains

!------------------------------------------------------------------

subroutine update_inflation(old_inflate, old_inflate_sd, new_inflate, new_inflate_sd, prior_mean, prior_var, &
   obs, obs_var, gamma)

! Given information of the old inflate and inflate_sd,
! the ensemble prior_mean and prior_var for an observation, and the observed value
! and observational error variance, computes updated values for the inflate and
! inflate_sd values using adaptive inflation algorithm 
! The gamma paramter gives the localized prior correlation times the localization

real(r_single),                    intent(in)    :: prior_mean, prior_var ! ensemble mean/var for an obs
real(r_single),                    intent(in)    :: obs, obs_var, gamma ! observed obs value and obs error
real(r_single),                    intent(inout) :: new_inflate, new_inflate_sd ! new inflation mean and sd
real(r_single),                    intent(in) :: old_inflate, old_inflate_sd ! old inflation mean and sd

new_inflate = old_inflate
new_inflate_sd = old_inflate_sd

! If the inflate_sd not positive, keep everything the same
if(old_inflate_sd <= 0.0_r_single) return

! Use bayes theorem to update
call bayes_cov_inflate(prior_mean, prior_var, obs, obs_var, old_inflate, &
   old_inflate_sd, gamma, new_inflate, new_inflate_sd, prior_inf_sd_lower_bound)

if(isnan(new_inflate))        new_inflate = old_inflate
if(isnan(new_inflate_sd))        new_inflate_sd = old_inflate_sd

! Make sure inflate satisfies constraints
if(new_inflate < prior_inf_lower_bound) new_inflate = prior_inf_lower_bound
if(new_inflate > prior_inf_upper_bound) new_inflate = prior_inf_upper_bound

! Make sure sd satisfies constraints
if(new_inflate_sd < prior_inf_sd_lower_bound) new_inflate_sd = prior_inf_sd_lower_bound

end subroutine update_inflation

!------------------------------------------------------------------

subroutine bayes_cov_inflate(x_p, sigma_p_2, y_o, sigma_o_2, lambda_mean, lambda_sd, &
   gamma, new_cov_inflate, new_cov_inflate_sd, sd_lower_bound_in)

real(r_single), intent(in)  :: x_p, sigma_p_2, y_o, sigma_o_2, lambda_mean, lambda_sd, gamma
real(r_single), intent(in)  :: sd_lower_bound_in
real(r_single), intent(out) :: new_cov_inflate, new_cov_inflate_sd

integer  :: i, mlambda_index(1)

real(r_single) :: new_1_sd, new_max, ratio, lambda_sd_2
real(r_single) :: dist_2, b, c, d, Q, R, disc, alpha, beta, cube_root_alpha, cube_root_beta, x
real(r_single) :: rrr, cube_root_rrr, angle, mx(3), sep(3), mlambda(3)


! If gamma is 0, nothing happens
if(gamma <= 0.0_r_single) then
   new_cov_inflate = lambda_mean
   new_cov_inflate_sd = lambda_sd
   return
endif

! Computation saver
lambda_sd_2 = lambda_sd**2
dist_2 = (x_p - y_o)**2

! Use ONLY the linear approximation, cubic solution below can be numerically
! unstable for extreme cases. Should look at this later.
!!!if(gamma > 0.99_r8) then
if(gamma > 1.01_r_single) then
  print *,'in adaptive_inflation_mod: gamma =',gamma

! FIXME: rectify the comment in the following line (about gamma
! and 1.0) and the test above, which clearly used to be true for 
! 1.0 but is not anymore.

! The solution of the cubic below only works if gamma is 1.0
! Can analytically find the maximum of the product: d/dlambda is a
! cubic polynomial in lambda**2; solve using cubic formula for real root
! Can write so that coefficient of x**3 is 1, other coefficients are:
   b = -1.0_r_single * (sigma_o_2 + sigma_p_2 * lambda_mean)
   c = lambda_sd_2 * sigma_p_2**2 / 2.0_r_single
   d = -1.0_r_single * (lambda_sd_2 * sigma_p_2**2 * dist_2) / 2.0_r_single

   Q = c - b**2 / 3
   R = d + (2 * b**3) / 27 - (b * c) / 3

   ! Compute discriminant, if this is negative have 3 real roots, else 1 real root
   disc = R**2 / 4 + Q**3 / 27

   if(disc < 0.0_r_single) then
      rrr = sqrt(-1.0 * Q**3 / 27)
      ! Note that rrr is positive so no problem for cube root
      cube_root_rrr = rrr ** (1.0 / 3.0)
      angle = acos(-0.5 * R / rrr)
      do i = 0, 2
         mx(i+1) = 2.0_r_single * cube_root_rrr * cos((angle + i * 2.0_r_single * pi) / 3.0_r_single) - b / 3.0_r_single
         mlambda(i + 1) = (mx(i + 1) - sigma_o_2) / sigma_p_2
         sep(i+1) = abs(mlambda(i + 1) - lambda_mean)
      end do
      ! Root closest to initial peak is appropriate
      mlambda_index = minloc(sep)
      new_cov_inflate = mlambda(mlambda_index(1))

   else
      ! Only one real root here, find it.

      ! Compute the two primary terms
      alpha = -R/2 + sqrt(disc)
      beta = R/2 + sqrt(disc)

      cube_root_alpha = abs(alpha) ** (1.0 / 3.0) * abs(alpha) / alpha
      cube_root_beta = abs(beta) ** (1.0 / 3.0) * abs(beta) / beta
   
      x = cube_root_alpha - cube_root_beta - b / 3.0

      ! This root is the value of x = theta**2
      new_cov_inflate = (x - sigma_o_2) / sigma_p_2

   endif

   ! Put in code to approximate the mode (new_cov_inflate)
else
   ! Approximate with Taylor series for likelihood term
   call linear_bayes(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd_2, gamma, &
      new_cov_inflate)
endif

! Bail out to save cost when lower bound is reached on lambda standard deviation
if(lambda_sd <= sd_lower_bound_in) then
   new_cov_inflate_sd = lambda_sd
else
   ! Compute by forcing a Gaussian fit at one positive SD
! First compute the new_max value for normalization purposes
   new_max = compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, &
      gamma, new_cov_inflate)

! Find value at a point one OLD sd above new mean value
   new_1_sd = compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma, &
      new_cov_inflate + lambda_sd)

   ! If either the numerator or denominator of the following computation 
   ! of 'ratio' is going to be zero (or almost so), return the original incoming
   ! inflation value.  The computation would have resulted in either Inf or NaN.
   if (abs(new_max) <= TINY(0.0_r_single) .or. abs(new_1_sd) <= TINY(0.0_r_single)) then
      new_cov_inflate_sd = lambda_sd
      return
   endif

   ratio = new_1_sd / new_max 

   ! Another error for numerical issues; if ratio is larger than 0.99, bail out
   if(ratio > 0.99_r_single) then
      new_cov_inflate_sd = lambda_sd
      return
   endif

   ! Can now compute the standard deviation consistent with this as
      ! sigma = sqrt(-x^2 / (2 ln(r))  where r is ratio and x is lambda_sd (distance from mean)
   new_cov_inflate_sd = sqrt( -1.0_r_single * lambda_sd_2 / (2.0_r_single * log(ratio)))

   ! Prevent an increase in the sd of lambda???
   ! For now, this is mostly countering numerical errors in this computation
   if(new_cov_inflate_sd > lambda_sd) new_cov_inflate_sd = lambda_sd

endif

end subroutine bayes_cov_inflate

!------------------------------------------------------------------

function compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma, lambda)

! Used to update density by taking approximate gaussian product

real(r_single)             :: compute_new_density
real(r_single), intent(in) :: dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma, lambda

real(r_single) :: theta_2, theta
real(r_single) :: exponent_prior, exponent_likelihood


! Compute probability of this lambda being correct
exponent_prior = (lambda - lambda_mean)**2 / (-2.0_r_single * lambda_sd**2)

! Compute probability that observation would have been observed given this lambda
theta_2 = (1.0_r_single + gamma * (sqrt(lambda) - 1.0_r_single))**2 * sigma_p_2 + sigma_o_2
theta = sqrt(theta_2)

exponent_likelihood = dist_2 / ( -2.0_r_single * theta_2)

! Compute the updated probability density for lambda
! Have 1 / sqrt(2 pi) twice, so product is 1 / (2 pi)
compute_new_density = exp(exponent_likelihood + exponent_prior) / &
   (2.0_r_single * pi * lambda_sd * theta)

end function compute_new_density


!---------------------------------------------------------------------

subroutine linear_bayes(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd_2, gamma, &
   new_cov_inflate)

real(r_single), intent(in)    :: dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd_2
real(r_single), intent(in)    :: gamma
real(r_single), intent(inout) :: new_cov_inflate

real(r_single) :: theta_bar_2, u_bar, like_exp_bar, v_bar, like_bar, like_prime, theta_bar
real(r_single) :: a, b, c, plus_root, minus_root, dtheta_dlambda
   
! Compute value of theta at current lambda_mean
theta_bar_2 = (1.0_r_single + gamma * (sqrt(lambda_mean) - 1.0_r_single))**2 * sigma_p_2 + sigma_o_2
theta_bar = sqrt(theta_bar_2)
! Compute constant coefficient for likelihood at lambda_bar
u_bar = 1.0_r_single / (sqrt(2.0_r_single * pi) * theta_bar)
! Compute exponent of likelihood at lambda_bar
like_exp_bar = dist_2 / (-2.0_r_single * theta_bar_2)
! Compute exponential part of likelihood at lambda_bar
v_bar = exp(like_exp_bar)
! Compute value of likelihood at current lambda_bar value
like_bar = u_bar * v_bar

! If like_bar goes to 0, can't do anything, so just keep current values
if(like_bar <= 0.0_r_single) then
   new_cov_inflate = lambda_mean
   return
endif

! Next compute derivative of likelihood at this point

! First compute d/dlambda of theta evaluated at lambda_mean
! Verified correct by finite difference, 1 January, 2006
dtheta_dlambda = 0.5_r_single * sigma_p_2 * gamma *(1.0_r_single - gamma + gamma*sqrt(lambda_mean)) / &
   (theta_bar * sqrt(lambda_mean))
like_prime = (u_bar * v_bar * dtheta_dlambda / theta_bar) * (dist_2 / theta_bar_2 - 1.0_r_single)

! If like_prime goes to 0, can't do anything, so just keep current values
if(like_prime == 0.0_r_single) then
   new_cov_inflate = lambda_mean
   return
endif

a = 1.0_r_single
b = like_bar / like_prime - 2.0_r_single * lambda_mean
c = lambda_mean**2 -lambda_sd_2 - like_bar * lambda_mean / like_prime

! Use nice scaled quadratic solver to avoid precision issues
call solve_quadratic(a, b, c, plus_root, minus_root)

! Do a check to pick closest root
if(abs(minus_root - lambda_mean) < abs(plus_root - lambda_mean)) then
   new_cov_inflate = minus_root
else
   new_cov_inflate = plus_root
endif

end subroutine linear_bayes


!------------------------------------------------------------------------

subroutine comp_likelihood(dist_2, sigma_p_2, sigma_o_2, lambda_mean, gamma, like_bar)
real(r_single), intent(in)  :: dist_2, sigma_p_2, sigma_o_2, lambda_mean, gamma
real(r_single), intent(out) :: like_bar

real :: theta_bar_2, u_bar, like_exp_bar, v_bar

! Compute value of theta at current lambda_bar
theta_bar_2 = (1.0_r_single + gamma * (sqrt(lambda_mean) - 1.0_r_single))**2 * sigma_p_2 + sigma_o_2
! Compute constanc coefficient for likelihood at lambda_bar
u_bar = 1.0_r_single / (sqrt(2.0_r_single * pi) * sqrt(theta_bar_2))
! Compute exponent of likelihood at lambda_bar
like_exp_bar = dist_2 / (-2.0_r_single * theta_bar_2)
! Compute exponential part of likelihood at lambda_bar
v_bar = exp(like_exp_bar)
! Compute value of likelihood at current lambda_bar value
like_bar = u_bar * v_bar

end subroutine comp_likelihood

!------------------------------------------------------------------------

subroutine solve_quadratic(a, b, c, r1, r2)

real(r_single), intent(in)  :: a, b, c
real(r_single), intent(out) :: r1, r2

real(r_single) :: scaling, as, bs, cs, disc

! Scale the coefficients to get better round-off tolerance
scaling = max(abs(a), abs(b), abs(c))
as = a / scaling
bs = b / scaling
cs = c / scaling

! Get discriminant of scaled equation
disc = sqrt(bs**2 - 4.0_r_single * as * cs)

if(bs > 0.0_r_single) then
   r1 = (-bs - disc) / (2 * as)
else
   r1 = (-bs + disc) / (2 * as)
endif

! Compute the second root given the larger one
r2 = (cs / as) / r1

end subroutine solve_quadratic

end module adaptive_inflate_mod
