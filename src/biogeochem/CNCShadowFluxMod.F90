module CNCShadowFluxMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for shadow carbon flux variable update, non-mortality fluxes.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use clm_varpar                         , only : ndecomp_cascade_transitions, nlevdecomp
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: CShadowFlux1

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CShadowFlux1(num_soilc, filter_soilc, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
       shadow_soilbiogeochem_carbonflux_inst, shadow_soilbiogeochem_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon shadow flux
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                              , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                              , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilbiogeochem_carbonflux_type) , intent(in)    :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type), intent(in)    :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: shadow_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type), intent(in)    :: shadow_soilbiogeochem_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: l,fc,cc,j
    integer :: cdp
    !-----------------------------------------------------------------------

    associate(                                                               &
         cascade_donor_pool       => decomp_cascade_con%cascade_donor_pool , &
         soilbiogeochem_cs        => soilbiogeochem_carbonstate_inst       , &
         soilbiogeochem_cf        => soilbiogeochem_carbonflux_inst        , &
         shadow_soilbiogeochem_cs => shadow_soilbiogeochem_carbonstate_inst, &
         shadow_soilbiogeochem_cf => shadow_soilbiogeochem_carbonflux_inst   &
         )

      ! column-level non-mortality fluxes

      do fc = 1, num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_cascade_transitions
               cdp = cascade_donor_pool(l)
               if ( soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) /= 0._r8) then
                  shadow_soilbiogeochem_cf%decomp_cascade_hr_vr_col(cc,j,l) = &
                      soilbiogeochem_cf%decomp_cascade_hr_vr_col(cc,j,l) * &
                      (shadow_soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) &
                         / soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp))
               else
                  shadow_soilbiogeochem_cf%decomp_cascade_hr_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

      do fc = 1, num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_cascade_transitions
               cdp = cascade_donor_pool(l)
               if ( soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) /= 0._r8) then
                  shadow_soilbiogeochem_cf%decomp_cascade_ctransfer_vr_col(cc,j,l) = &
                      soilbiogeochem_cf%decomp_cascade_ctransfer_vr_col(cc,j,l) * &
                      (shadow_soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) &
                         / soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp))
               else
                  shadow_soilbiogeochem_cf%decomp_cascade_ctransfer_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

    end associate

  end subroutine CShadowFlux1

end module CNCShadowFluxMod
