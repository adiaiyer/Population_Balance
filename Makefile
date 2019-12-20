# This is the main makefile.
# Use with GNU make.
# Relies on makedepf90 to do dependency lists:
# You can obtain this free program from:
#   http://www.helsinki.fi/~eedelman/makedepf90.html
# New $ARCH entries may be required, depending on your compiler.
# Also requires fpx3 fortran preprocessor available at
#   http://wwwuser.gwdg.de/~jbehren/fpx3.html

MK_INCL_PATH=.
include Makefile.in

EXE := lesgo

SRCS =  bottombc.f90 \
        convec_crossflow.f90 \
        ddx.f90 \
        ddxy.f90 \
        ddy.f90 \
        ddz_uv.f90 \
        ddz_w.f90 \
        dealias1.f90 \
        dealias2.f90 \
        debug_mod.f90 \
        divstress_uv.f90 \
        divstress_w.f90 \
        dns_stress.f90 \
        energy.f90 \
        fft.f90 \
        filt_da.f90 \
        forcing.f90 \
        ic.f90 \
        ic_dns.f90 \
        immersedbc.f90 \
        initial_multi_crossflow_fringe.f90 \
        interpolag_Sdep_crossflow.f90 \
        interpolag_Ssim.f90 \
        io_crossflow_mod.f90 \
        lagrange_Sdep.f90 \
        lagrange_Ssim.f90 \
        main_crossflow_Coriol3D.f90 \
        messages.f90 \
        padd.f90 param.f90 \
        press_stag_array.f90 \
        ran3.f90 rmsdiv.f90 \
        scaledep_dynamic.f90 \
        scalars_module_crossflow_mod.f90 \
        scalars_module2_lab.f90 \
        sgs_stag.f90 \
        sgsmodule.f90 \
        sim_param.f90 \
        std_dynamic.f90 \
        string_util.f90 \
        test_filtermodule.f90 \
        topbc.f90 \
        tridag.f90 \
        tridag_array.f90 \
        types.f90 \
        unpadd.f90 \
        wallstress.f90 \
        wallstress_dns.f90 \
        spectral_interp.f90 \
        box_filter.f90 \
        dynamic_sc.f90 \
        outflow_fringe.f90

#SRCS =  sim_param.f90 \
	param_crossflow_Coriol3D.f90 \
        types.f90 \
        main_crossflow_Coriol3D.f90 \
	bottombc.f90 \
        convec_crossflow.f90 \
        ddx.f90 \
        ddxy.f90 \
        ddy.f90 \
        ddz_uv.f90 \
        ddz_w.f90 \
        dealias1.f90 \
        dealias2.f90 \
        debug_mod.f90 \
        divstress_uv.f90 \
        divstress_w.f90 \
        dns_stress.f90 \
        energy.f90 \
        fft.f90 \
        filt_da.f90 \
        forcing.f90 \
        ic.f90 \
        ic_dns.f90 \
        immersedbc.f90 \
        initial_multi.f90 \
        interpolag_Sdep_crossflow.f90 \
        interpolag_Ssim.f90 \
        io_crossflow.f90 \
        lagrange_Sdep.f90 \
        lagrange_Ssim.f90 \
        messages.f90 \
        padd.f90 \
        press_stag_array.f90 \
        ran3.f90 \
	rmsdiv.f90 \
        scaledep_dynamic.f90 \
        scalars_module_crossflow.f90 \
        scalars_module2.f90 \
        sgs_stag.f90 \
        sgsmodule.f90 \
        std_dynamic.f90 \
        string_util.f90 \
        test_filtermodule.f90 \
        topbc.f90 \
        tridag.f90 \
        tridag_array.f90 \
        unpadd.f90 \
        wallstress.f90 \
        wallstress_dns.f90 \
        spectral_interp.f90 \
        box_filter.f90 \
        dynamic_sc.f90

#SRCS =  sim_param.f90 \
        param_crossflow_Coriol3D.f90 \
        types.f90 \
        main_crossflow_Coriol3D.f90 \
        io_crossflow.f90 \
        sgsmodule.f90 \
        scalars_module_crossflow.f90 \
        bottombc.f90 \
        test_filtermodule.f90 \
        fft.f90 \
        topbc.f90 \
        scalars_module2.f90 \
        dealias1.f90 \
        dealias2.f90 \
        ddx.f90 \
        ddxy.f90 \
        ddy.f90 \
        ddz_uv.f90 \
        ddz_w.f90 \
        ran3.f90 \
        padd.f90 \
        unpadd.f90 \
        filt_da.f90 \
        immersedbc.f90 \
        debug_mod.f90 \
        initial_multi.f90 \
        wallstress.f90 \
        wallstress_dns.f90 \
        sgs_stag.f90 \
        divstress_uv.f90 \
        divstress_w.f90 \
        convec_crossflow.f90 \
        press_stag_array.f90 \
        forcing.f90 \
        energy.f90 \
        rmsdiv.f90 \
        lagrange_Sdep.f90 \
        lagrange_Ssim.f90 \
        interpolag_Sdep_crossflow.f90 \
        interpolag_Ssim.f90 \


LVLSET_SRCS = level_set_base.f90 level_set.f90 linear_simple.f90

CYL_SKEW_LS_SRCS = cyl_skew_base_ls.f90 cyl_skew_ls.f90

RNS_LS_SRCS = rns_base_ls.f90 rns_ls.f90

RNS_CYL_SKEW_LS_SRCS = rns_cyl_skew_ls.f90

TURBINES_SRCS = turbines.f90 turbines_base.f90

CPS_SRCS = concurrent_precursor.f90

ifeq ($(USE_MPI), yes)
  SRCS += mpi_transpose_mod.f90 tridag_array_pipelined.f90 #mpi_defs.f90
  EXE := $(EXE)-mpi
endif

ifeq ($(USE_CPS), yes)
  SRCS += $(CPS_SRCS)
  EXE := $(EXE)-cps
endif

ifeq ($(USE_LVLSET), yes)
  SRCS += $(LVLSET_SRCS)
  EXE := $(EXE)-ls
endif

ifeq ($(USE_RNS_LS), yes)
  SRCS += $(RNS_LS_SRCS)
  ifeq ($(USE_CYL_SKEW_LS), yes)
    SRCS += $(RNS_CYL_SKEW_LS_SRCS)
  endif
  EXE := $(EXE)-rns
endif

ifeq ($(USE_CYL_SKEW_LS), yes)
  SRCS += $(CYL_SKEW_LS_SRCS)
  EXE := $(EXE)-cs
endif

ifeq ($(USE_TURBINES), yes)
  SRCS += $(TURBINES_SRCS)
  EXE := $(EXE)-turbines
endif

ifeq ($(OUTPUT_EXTRA), yes)
  EXE := $(EXE)-exout
endif

ifeq ($(USE_DYN_TN), yes)
  EXE := $(EXE)-dyntn
endif

ifeq ($(USE_BINARY), yes)
  EXE := $(EXE)-binary
endif

ifeq ($(USE_SAFETYMODE), no)
  EXE := $(EXE)-safety_off
endif


#COMPSTR = '$(FPP) $$< > t.$$<; $$(FC) -c -o $$@ $$(FFLAGS) t.$$<; rm -f t.$$<'
COMPSTR = '$(FPP) $$< > t.$$<; $$(FC) -c -o $$@ $$(FFLAGS) t.$$<'


include .depend

.depend: $(SRCS)
	mkdir -p $(OPATH) $(MPATH);
	makedepf90 -r $(COMPSTR) -b $(OPATH) -o $(EXE) $(SRCS) > .depend

debug:
	$(MAKE) $(EXE) "FFLAGS = $(FDEBUG) $(FFLAGS)"

prof:
	$(MAKE) $(EXE) "FFLAGS = $(FPROF) $(FFLAGS)"

# Other support programs are listed below this point
convert_endian:	utils/convert_endian.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

# This doesn't remove .mod files--should be OK as long a dependency list 
# for the .o files is correct.
# FOBJ is defined in .depend
.PHONY : clean
clean :
	rm -rf $(OPATH)/* $(FOBJ) .depend* $(MPATH)/*.mod
	rm -f t.*
