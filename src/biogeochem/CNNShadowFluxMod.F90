module CNNShadowFluxMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for shadow nitrogen flux variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod                         , only : r8 => shr_kind_r8
  use clm_varpar                           , only : ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools
  use SoilBiogeochemCarbonFluxType         , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType      , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType       , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemDecompCascadeConType   , only : decomp_cascade_con
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: NShadowFlux1
  public  :: NShadowFlux3

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NShadowFlux1(num_soilc, filter_soilc, shadow_soilbiogeochem_carbonflux_inst, &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst, &
       shadow_soilbiogeochem_nitrogenflux_inst, shadow_soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the nitrogen shadow flux
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilbiogeochem_carbonflux_type)   , intent(in)    :: shadow_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(in)    :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: shadow_soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: shadow_soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: l,fc,cc,j
    integer  :: cdp, crp
    real(r8) :: shadow_base_ratio
    !-----------------------------------------------------------------------

    associate(                                                                                &
         cascade_donor_pool             => decomp_cascade_con%cascade_donor_pool            , &
         cascade_receiver_pool          => decomp_cascade_con%cascade_receiver_pool         , &
         floating_cn_ratio_decomp_pools => decomp_cascade_con%floating_cn_ratio_decomp_pools, &
         initial_cn_ratio               => decomp_cascade_con%initial_cn_ratio              , &
         soilbiogeochem_ns              => soilbiogeochem_nitrogenstate_inst                , &
         soilbiogeochem_nf              => soilbiogeochem_nitrogenflux_inst                 , &
         shadow_soilbiogeochem_cf       => shadow_soilbiogeochem_carbonflux_inst            , &
         shadow_soilbiogeochem_ns       => shadow_soilbiogeochem_nitrogenstate_inst         , &
         shadow_soilbiogeochem_nf       => shadow_soilbiogeochem_nitrogenflux_inst            &
         )

      ! column-level non-mortality fluxes

      do fc = 1, num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_cascade_transitions
               cdp = cascade_donor_pool(l)
               if (soilbiogeochem_ns%decomp_npools_vr_col(cc,j,cdp) /= 0._r8) then
                  shadow_base_ratio = shadow_soilbiogeochem_ns%decomp_npools_vr_col(cc,j,cdp) &
                     / soilbiogeochem_ns%decomp_npools_vr_col(cc,j,cdp)
                  shadow_soilbiogeochem_nf%decomp_cascade_ntransfer_vr_col(cc,j,l) = &
                     soilbiogeochem_nf%decomp_cascade_ntransfer_vr_col(cc,j,l) * shadow_base_ratio
                  crp = cascade_receiver_pool(l)
                  if (floating_cn_ratio_decomp_pools(crp)) then
                     ! set to zero for now
                     ! this is not applicable to MIMICS
                     shadow_soilbiogeochem_nf%decomp_cascade_sminn_flux_vr_col(cc,j,l) = 0._r8
                  else
                     shadow_soilbiogeochem_nf%decomp_cascade_sminn_flux_vr_col(cc,j,l) = &
                        shadow_soilbiogeochem_cf%decomp_cascade_ctransfer_vr_col(cc,j,l) &
                        * (1._r8 / initial_cn_ratio(crp)) &
                        - shadow_soilbiogeochem_nf%decomp_cascade_ntransfer_vr_col(cc,j,l)
                  end if
               else
                  shadow_soilbiogeochem_nf%decomp_cascade_ntransfer_vr_col(cc,j,l) = 0._r8
                  shadow_soilbiogeochem_nf%decomp_cascade_sminn_flux_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

    end associate

  end subroutine NShadowFlux1

  !-----------------------------------------------------------------------
  subroutine NShadowFlux3(num_soilc, filter_soilc, &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst, &
       shadow_soilbiogeochem_nitrogenflux_inst, shadow_soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon shadow fluxes for fire mortality
    !
    ! !ARGUMENTS:
    integer                                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilbiogeochem_nitrogenflux_type) , intent(in)    :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: shadow_soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: shadow_soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: l,fc,cc,j
    integer  :: cdp, crp
    real(r8) :: shadow_base_ratio
    !-----------------------------------------------------------------------

    associate(                                                                                &
         soilbiogeochem_ns              => soilbiogeochem_nitrogenstate_inst                , &
         soilbiogeochem_nf              => soilbiogeochem_nitrogenflux_inst                 , &
         shadow_soilbiogeochem_ns       => shadow_soilbiogeochem_nitrogenstate_inst         , &
         shadow_soilbiogeochem_nf       => shadow_soilbiogeochem_nitrogenflux_inst            &
         )

      do fc = 1, num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_pools
               if ( soilbiogeochem_ns%decomp_npools_vr_col(cc,j,l) /= 0._r8) then
                  shadow_base_ratio = shadow_soilbiogeochem_ns%decomp_npools_vr_col(cc,j,l) &
                     / soilbiogeochem_ns%decomp_npools_vr_col(cc,j,l)
                  shadow_soilbiogeochem_nf%m_decomp_npools_to_fire_vr_col(cc,j,l)  =  &
                     soilbiogeochem_nf%m_decomp_npools_to_fire_vr_col(cc,j,l) * shadow_base_ratio
               else
                  shadow_soilbiogeochem_nf%m_decomp_npools_to_fire_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

    end associate

  end subroutine NShadowFlux3

end module CNNShadowFluxMod
