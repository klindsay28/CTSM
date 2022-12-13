module CNNShadowFluxMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for shadow nitrogen flux variable update, non-mortality fluxes.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use clm_varpar                         , only : ndecomp_cascade_transitions, nlevdecomp
  use SoilBiogeochemNitrogenStateType      , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType       , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemDecompCascadeConType   , only : decomp_cascade_con
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: NShadowFlux1

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NShadowFlux1(num_soilc, filter_soilc, &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst, &
       shadow_soilbiogeochem_nitrogenflux_inst, shadow_soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the nitrogen shadow flux
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                              , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                              , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilbiogeochem_nitrogenflux_type) , intent(in)    :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: shadow_soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: shadow_soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: l,fc,cc,j
    integer :: cdp
    !-----------------------------------------------------------------------

    associate(                                                               &
         cascade_donor_pool       => decomp_cascade_con%cascade_donor_pool , &
         soilbiogeochem_ns        => soilbiogeochem_nitrogenstate_inst       , &
         soilbiogeochem_nf        => soilbiogeochem_nitrogenflux_inst        , &
         shadow_soilbiogeochem_ns => shadow_soilbiogeochem_nitrogenstate_inst, &
         shadow_soilbiogeochem_nf => shadow_soilbiogeochem_nitrogenflux_inst   &
         )

      ! column-level non-mortality fluxes

      do fc = 1, num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_cascade_transitions
               cdp = cascade_donor_pool(l)
               if ( soilbiogeochem_ns%decomp_npools_vr_col(cc,j,cdp) /= 0._r8) then
                  shadow_soilbiogeochem_nf%decomp_cascade_ntransfer_vr_col(cc,j,l) = &
                      soilbiogeochem_nf%decomp_cascade_ntransfer_vr_col(cc,j,l) * &
                      (shadow_soilbiogeochem_ns%decomp_npools_vr_col(cc,j,cdp) &
                         / soilbiogeochem_ns%decomp_npools_vr_col(cc,j,cdp))
               else
                  shadow_soilbiogeochem_nf%decomp_cascade_ntransfer_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

      do fc = 1, num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_cascade_transitions
               cdp = cascade_donor_pool(l)
               if ( soilbiogeochem_ns%decomp_npools_vr_col(cc,j,cdp) /= 0._r8) then
                  shadow_soilbiogeochem_nf%decomp_cascade_sminn_flux_vr_col(cc,j,l) = &
                      soilbiogeochem_nf%decomp_cascade_sminn_flux_vr_col(cc,j,l) * &
                      (shadow_soilbiogeochem_ns%decomp_npools_vr_col(cc,j,cdp) &
                         / soilbiogeochem_ns%decomp_npools_vr_col(cc,j,cdp))
               else
                  shadow_soilbiogeochem_nf%decomp_cascade_sminn_flux_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

    end associate

  end subroutine NShadowFlux1

end module CNNShadowFluxMod
