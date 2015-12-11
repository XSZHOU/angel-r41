# MAKEFILE PROJECT ANGEL
# independent of platform and compiler
#
# syntax: targetfile : sourcefiles
#                      commands

BASE_SRC_PATH=$(shell pwd)

all: make.inc umfpack flens binary #rappture # S.Z.

# all: make.inc library binary

# Added by Steve Clark, Nanohub
# only used by Nanohub team
install:
	. /etc/environ.sh ; use -e -r mpich2-1.3-gnu-4.7.2; make all


make.inc: set_make_inc.sh
	./set_make_inc.sh $(BASE_SRC_PATH)


include make.inc


    
OBJFILES = obj/all.o \
    obj/Exception.o              obj/Logger.o                 obj/Timer.o                 obj/Filenames.o        \
    obj/MPIs.o                   obj/Interrupt.o              obj/MatFlens.o              obj/MatBanded.o        \
    obj/Vertex.o                 obj/Edge.o                   obj/Face.o                  obj/Element.o          \
    obj/Region.o                 obj/Contact.o                obj/Geometry.o              obj/BoxMethod.o        \
    obj/DomainPoint.o            obj/DomainNode.o             obj/DomainMaster.o          obj/Domain.o           \
    obj/InputParser.o            obj/MaterialDatabase.o       obj/StrainPolarization.o    obj/Options.o          \
    obj/PoissonFEM1D.o           obj/PoissonProblem.o         obj/Voltages.o                                     \
    obj/CSRMatrix.o              obj/LinearSolverUmfpack.o    \
    obj/Kspace.o                 obj/Energies.o               obj/Hamiltonian.o           obj/TdkpInfoDesk.o     \
    obj/InterfaceEffMass.o       obj/InterfaceEffMassOrtho.o  obj/InterfaceEffMassOrtho2Band.o obj/Overlap.o     \
    obj/NEGFObject.o             obj/GreenFunctions.o         obj/ContactFermilevel.o                            \
    obj/SelfEnergy.o             obj/SEContacts.o             obj/SEBuettiker.o           obj/SEGolizadeh.o      \
    obj/SEOpticalPhonon.o        obj/SEAcousticPhonon.o       obj/SEIonizedImpurities.o                          \
    obj/SEPhotonSpontaneous.o    obj/Luminescence.o           obj/Transmission.o          obj/QuasiFermilevel.o  \
    obj/SelfEnergies.o           obj/Resonances.o             obj/InnerLoop.o             obj/Current.o          \
    obj/PostProcessing.o         obj/OuterLoop.o       


# files corresponding to includes in all.h
ALL_H    = common/all.h common/Logger.h common/Exception.h common/Filenames.h common/Timer.h common/Interrupt.h \
			common/MPIs.h common/Mat.h common/MatFlens.inl common/MatBanded.inl
# everything related to geometry
GEOM_H = geometry/Geometry.h geometry/Contact.h geometry/Region.h geometry/Element.h geometry/Face.h geometry/Edge.h geometry/Vertex.h
# k-space, E-space
KSPACE_ESPACE = kspace_Espace/Energies.h kspace_Espace/Kspace.h

	
clean:
	rm -f obj/*.o negf_wrap.cpp
	rm -f ../rappture/*.o
	rm -f angel.bin
	
distclean: clean
	rm -f ../lib/*.so
	(cd umfpack/UFconfig; rm -f UFconfig.mk; ln -s $(UF_MK) UFconfig.mk)
	(cd umfpack/AMD; make clean)
	(cd umfpack/UMFPACK; make clean)
	rm -f UFconfig.mk
	(cd flens; make clean)
	rm -rf flens/flens/.obj
	rm -rf flens/flens/.dep
	rm -f make.inc

umfpack: ../lib/libumfpack.so


../lib/libumfpack.so: umfpack/*
	@echo "Compiling UMFPACK..."
	(cd umfpack/UFconfig; rm -f UFconfig.mk; ln -s $(UF_MK) UFconfig.mk)
	(cd umfpack/AMD; make library)
	(cd umfpack/UMFPACK; make library)
	(cd umfpack; ./build_shared)
	cp umfpack/libumfpack.so ../lib/

flens: ../lib/libflens.so
# notice here, change the CXX to be g++ in /flens/flens/Makefile, otherwise linking error
../lib/libflens.so: flens/* 
	@echo "Compiling FLENS..."
	(cd flens/flens; make FLENS_HOME="$(PWD)/flens")
	cp flens/libflens.so ../lib/

binary: $(OBJFILES) $(OMPDUMMY) obj/main.o
	$(CPP) $(CPP_FLAGS) $(OBJFILES) $(OMPDUMMY) obj/main.o -o ../bin/angel.bin $(CPP_LIBS)


#S.Z. ok! now the part below is not relevant...

ifeq (1,1)

rappture: ../rappture/rapp_angel.bin ../rappture/see_results.bin


../rappture/rapp_angel.bin: ../rappture/rapp_angel.cpp
	(cd ../rappture; \
             $(CPP) $(RAPP_INC) -DNEGF_DIR=\"$(NEGF_VAR)\" -c rapp_angel.cpp; \
             $(CPP) -o rapp_angel.bin rapp_angel.o $(RAPP_LIB))

../rappture/see_results.bin: ../rappture/see_results.cpp
	(cd ../rappture; \
             $(CPP) $(RAPP_INC) -c see_results.cpp -DNEGF_DIR=\"$(NEGF_VAR)\"; \
             $(CPP) -o see_results.bin see_results.o $(RAPP_LIB))


library: $(OBJFILES) obj/negf_wrap.o
	$(CPP) $(CPP_FLAGS) -shared $(OBJFILES) obj/negf_wrap.o -o _negf.so $(CPP_LIBS)

mattest: obj/all.o obj/Exception.o obj/Logger.o obj/Timer.o obj/Filenames.o obj/MPIs.o obj/Interrupt.o obj/MatFlens.o obj/MatBanded.o obj/mattest.o
	$(CPP) $(CPP_FLAGS) obj/all.o obj/Exception.o obj/Logger.o obj/Timer.o obj/Filenames.o obj/MPIs.o obj/Interrupt.o obj/MatFlens.o obj/MatBanded.o \
	obj/mattest.o -o mattest.bin $(CPP_LIBS)
	
commtest: obj/all.o obj/Exception.o obj/Logger.o obj/Timer.o obj/Filenames.o obj/MPIs.o obj/Interrupt.o obj/MatFlens.o obj/MatBanded.o obj/commtest.o
	$(CPP) $(CPP_FLAGS) obj/all.o obj/Exception.o obj/Logger.o obj/Timer.o obj/Filenames.o obj/MPIs.o obj/Interrupt.o obj/MatFlens.o obj/MatBanded.o \
	obj/commtest.o -o commtest.bin $(CPP_LIBS)

endif

#=====================================================================================================	

obj/negf_wrap.o: common/*.h negf/*.h geometry/*.h kspace_Espace/*.h io/*.h newton/*.h poisson/*.h \
			 negf.i negf_wrap.cpp
	$(CPP) -Wall -fPIC -g -m64 -O -c $(DEFINES) negf_wrap.cpp -o obj/negf_wrap.o -I$(PYTHON_H) $(CPP_ALL_INC)

negf_wrap.cpp: common/*.h negf/*.h geometry/*.h kspace_Espace/*.h io/*.h newton/*.h poisson/*.h \
			   negf.i
	$(SWIG) $(DEFINES) $(CPP_INC) -I$(PYTHON_H) -c++ -python -o negf_wrap.cpp negf.i

obj/main.o: */*.h */*/*.h main.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/main.o main.cpp
	
obj/ompdummy.o: ompdummy.c
	cc -c $(CPP_FLAGS) -o obj/ompdummy.o ompdummy.c

obj/mattest.o: */*.h */*/*.h common/*.inl mattest.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/mattest.o mattest.cpp

obj/commtest.o: */*.h */*/*.h common/*.inl commtest.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/commtest.o commtest.cpp

obj/all.o: $(ALL_H) common/all.cpp  
	$(CPP) -c $(CPP_FLAGS) -DNEGF_DIR=\"$(NEGF_VAR)\" $(CPP_ALL_INC) -o obj/all.o common/all.cpp
		
obj/Exception.o: common/Exception.cpp common/Exception.h  
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Exception.o common/Exception.cpp
	
obj/Logger.o: $(ALL_H) common/Logger.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Logger.o common/Logger.cpp
	
obj/Timer.o: $(ALL_H) common/Timer.cpp 
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Timer.o common/Timer.cpp
	
obj/Filenames.o: $(ALL_H) common/Filenames.cpp 
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Filenames.o common/Filenames.cpp
	
obj/MPIs.o: $(ALL_H) common/MPIs.cpp 
	$(CPP) -c $(CPP_FLAGS) $(CPP_ALL_INC) -o obj/MPIs.o common/MPIs.cpp
	
obj/Interrupt.o: $(ALL_H) common/Interrupt.cpp 
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Interrupt.o common/Interrupt.cpp
	
obj/MatFlens.o: $(ALL_H) common/MatFlens.cpp common/MatFlens.inl common/MatBanded.inl
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/MatFlens.o common/MatFlens.cpp
	
obj/MatBanded.o: $(ALL_H) common/MatBanded.cpp common/MatFlens.inl common/MatBanded.inl
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/MatBanded.o common/MatBanded.cpp
	
obj/Vertex.o: $(ALL_H) geometry/Contact.h \
			  geometry/Vertex.cpp geometry/Vertex.h         
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Vertex.o geometry/Vertex.cpp
	
obj/Edge.o: $(ALL_H) geometry/Vertex.h \
			geometry/Edge.cpp geometry/Edge.h
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Edge.o geometry/Edge.cpp
	
obj/Face.o: $(ALL_H) geometry/Edge.h geometry/Vertex.h \
			geometry/Face.h geometry/Face.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Face.o geometry/Face.cpp

obj/Element.o: $(ALL_H) geometry/Vertex.h geometry/Edge.h geometry/Face.h \
			   geometry/Element.cpp geometry/Element.h
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Element.o geometry/Element.cpp

obj/Region.o: $(ALL_H) io/PropertyContainer.h geometry/Element.h \
			  geometry/Region.h geometry/Region.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Region.o geometry/Region.cpp
	
obj/Contact.o: $(ALL_H) geometry/Vertex.h geometry/Region.h \
			   geometry/Contact.h geometry/Contact.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Contact.o geometry/Contact.cpp
	
obj/Geometry.o: $(ALL_H) geometry/Vertex.h geometry/Edge.h geometry/Element.h \
				geometry/Region.h geometry/Face.h \
				geometry/Contact.h geometry/Geometry.h geometry/Geometry.cpp 
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Geometry.o geometry/Geometry.cpp   
	 
obj/InputParser.o: $(ALL_H) $(GEOM_H) io/InputParser.h io/InputParser.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_ALL_INC) -o obj/InputParser.o io/InputParser.cpp
	
obj/OutputData.o: $(ALL_H) $(GEOM_H) io/InputParser.h newton/Equation.h \
				  io/OutputData.h io/OutputData.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/OutputData.o io/OutputData.cpp
	
obj/BoxMethod.o: $(ALL_H) $(GEOM_H) geometry/BoxMethod.h geometry/BoxMethod.cpp 
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/BoxMethod.o geometry/BoxMethod.cpp    
	
obj/Overlap.o: $(ALL_H) $(GEOM_H) io/Options.h \
				geometry/Overlap.h geometry/Overlap.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Overlap.o geometry/Overlap.cpp 
	
obj/DomainPoint.o: $(ALL_H) kspace_Espace/DomainPoint.h kspace_Espace/DomainPoint.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/DomainPoint.o kspace_Espace/DomainPoint.cpp 

obj/DomainNode.o: $(ALL_H) kspace_Espace/DomainPoint.h \
				kspace_Espace/DomainNode.h kspace_Espace/DomainNode.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/DomainNode.o kspace_Espace/DomainNode.cpp 

obj/DomainMaster.o: $(ALL_H) kspace_Espace/DomainNode.h \
				kspace_Espace/DomainMaster.h kspace_Espace/DomainMaster.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/DomainMaster.o kspace_Espace/DomainMaster.cpp 

obj/Domain.o: $(ALL_H) kspace_Espace/DomainMaster.h \
				kspace_Espace/Domain.h kspace_Espace/Domain.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Domain.o kspace_Espace/Domain.cpp 
	
obj/Kspace.o: $(ALL_H) kspace_Espace/DomainMaster.h kspace_Espace/Domain.h \
				io/Options.h kspace_Espace/DomainPoint.h \
				io/PropertyContainer.h io/TernaryPropertyContainer.h io/MaterialDatabase.h \
				kspace_Espace/Kspace.h kspace_Espace/Kspace.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Kspace.o kspace_Espace/Kspace.cpp 
	
obj/Energies.o: $(ALL_H) kspace_Espace/DomainMaster.h kspace_Espace/Domain.h \
				io/Options.h kspace_Espace/DomainPoint.h \
				kspace_Espace/Energies.h kspace_Espace/Energies.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Energies.o kspace_Espace/Energies.cpp 
	
obj/Hamiltonian.o: $(ALL_H) $(GEOM_H) kspace_Espace/DomainMaster.h kspace_Espace/Domain.h \
				io/Options.h kspace_Espace/DomainPoint.h io/PropertyContainer.h \
				io/TernaryPropertyContainer.h io/MaterialDatabase.h poisson/StrainPolarization.h \
				negf/hamiltonian/*.h negf/hamiltonian/Hamiltonian.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/Hamiltonian.o negf/hamiltonian/Hamiltonian.cpp 
	
obj/ContactFermilevel.o: $(ALL_H) $(GEOM_H) kspace_Espace/DomainMaster.h kspace_Espace/Domain.h \
				io/Options.h kspace_Espace/DomainPoint.h io/PropertyContainer.h \
				negf/hamiltonian/Hamiltonian.h negf/selfenergies/SelfEnergy.h negf/NEGFObject.h \
				negf/selfenergies/SEContacts.h negf/ContactFermilevel.h negf/ContactFermilevel.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o obj/ContactFermilevel.o negf/ContactFermilevel.cpp 
	
obj/TdkpInfoDesk.o: $(ALL_H) $(GEOM_H) io/PropertyContainer.h io/TernaryPropertyContainer.h\
				io/MaterialDatabase.h poisson/StrainPolarization.h \
				negf/hamiltonian/TdkpInfoDesk.h negf/hamiltonian/TdkpInfoDesk.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) $(TDKP_INC) -o obj/TdkpInfoDesk.o negf/hamiltonian/TdkpInfoDesk.cpp 
	
obj/InterfaceEffMass.o: $(ALL_H) $(GEOM_H) negf/hamiltonian/InterfaceNEGFWell.h poisson/StrainPolarization.h \
				negf/hamiltonian/InterfaceEffMass.h negf/hamiltonian/InterfaceEffMass.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) $(TDKP_INC) -o obj/InterfaceEffMass.o negf/hamiltonian/InterfaceEffMass.cpp 
	
obj/InterfaceEffMassOrtho.o: $(ALL_H) $(GEOM_H) negf/hamiltonian/InterfaceNEGFWell.h poisson/StrainPolarization.h \
				negf/hamiltonian/InterfaceEffMassOrtho.h negf/hamiltonian/InterfaceEffMassOrtho.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) $(TDKP_INC) -o obj/InterfaceEffMassOrtho.o negf/hamiltonian/InterfaceEffMassOrtho.cpp 
	
obj/InterfaceEffMassOrtho2Band.o: $(ALL_H) $(GEOM_H) negf/hamiltonian/InterfaceNEGFWell.h poisson/StrainPolarization.h \
				negf/hamiltonian/InterfaceEffMassOrtho2Band.h negf/hamiltonian/InterfaceEffMassOrtho2Band.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) $(TDKP_INC) -o obj/InterfaceEffMassOrtho2Band.o negf/hamiltonian/InterfaceEffMassOrtho2Band.cpp 

# PropertyContainer        is template!
# TernaryPropertyContainer is template!	

obj/MaterialDatabase.o: $(ALL_H) io/PropertyContainer.h io/TernaryPropertyContainer.h\
						io/MaterialDatabase.h io/MaterialDatabase.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ io/MaterialDatabase.cpp  

obj/StrainPolarization.o: $(ALL_H) $(GEOM_H) io/PropertyContainer.h io/TernaryPropertyContainer.h\
						io/MaterialDatabase.h poisson/StrainPolarization.cpp poisson/StrainPolarization.h
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ poisson/StrainPolarization.cpp 

obj/Options.o: $(ALL_H) io/PropertyContainer.h io/TernaryPropertyContainer.h\
						io/InputParser.h io/Options.cpp io/Options.h
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ io/Options.cpp 
		
obj/PoissonFEM1D.o: $(ALL_H) $(GEOM_H) poisson/PoissonFEM1D.h poisson/PoissonFEM1D.cpp 
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ poisson/PoissonFEM1D.cpp	
	
obj/PoissonProblem.o: $(ALL_H) $(GEOM_H) io/MaterialDatabase.h \
				poisson/PoissonProblem.h poisson/PoissonProblem.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ poisson/PoissonProblem.cpp	

obj/CSRMatrix.o: $(ALL_H) poisson/CSRMatrix.h poisson/CSRMatrix.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ poisson/CSRMatrix.cpp
	
obj/LinearSolverUmfpack.o: $(ALL_H) poisson/CSRMatrix.h  \
					poisson/LinearSolverUmfpack.h poisson/LinearSolverUmfpack.cpp 
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) $(UMFPACK_INC) -o $@ poisson/LinearSolverUmfpack.cpp
		
obj/NEGFObject.o: $(ALL_H) $(GEOM_H) kspace_Espace/Energies.h kspace_Espace/Kspace.h \
                 io/Options.h negf/NEGFObject.h negf/NEGFObject.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/NEGFObject.cpp
		
obj/GreenFunctions.o: $(ALL_H) $(GEOM_H) kspace_Espace/Energies.h kspace_Espace/Kspace.h \
				 negf/NEGFObject.h geometry/Overlap.h \
                 negf/GreenFunctions.h negf/GreenFunctions.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/GreenFunctions.cpp
	
obj/SelfEnergy.o: $(ALL_H) $(GEOM_H) kspace_Espace/Energies.h kspace_Espace/Kspace.h \
				 negf/hamiltonian/Hamiltonian.h negf/NEGFObject.h \
                 negf/selfenergies/SelfEnergy.h negf/selfenergies/SelfEnergy.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/selfenergies/SelfEnergy.cpp
		
obj/SelfEnergies.o: $(ALL_H) $(GEOM_H) kspace_Espace/Energies.h kspace_Espace/Kspace.h \
				 geometry/Overlap.h negf/NEGFObject.h negf/selfenergies/*.h \
				 negf/hamiltonian/Hamiltonian.h negf/GreenFunctions.h \
                 negf/SelfEnergies.h negf/SelfEnergies.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/SelfEnergies.cpp
		
obj/SEContacts.o: $(ALL_H) $(GEOM_H) $(KSPACE_ESPACE) negf/NEGFObject.h negf/selfenergies/SelfEnergy.h \
				io/Options.h negf/hamiltonian/Hamiltonian.h geometry/Overlap.h \
                negf/selfenergies/SEContacts.h negf/selfenergies/SEContacts.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/selfenergies/SEContacts.cpp
		
obj/SEBuettiker.o: $(ALL_H) $(GEOM_H) $(KSPACE_ESPACE) negf/NEGFObject.h negf/selfenergies/SelfEnergy.h \
				io/Options.h negf/hamiltonian/Hamiltonian.h geometry/Overlap.h negf/GreenFunctions.h \
                negf/selfenergies/SEBuettiker.h negf/selfenergies/SEBuettiker.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/selfenergies/SEBuettiker.cpp
		
obj/SEGolizadeh.o: $(ALL_H) $(GEOM_H) $(KSPACE_ESPACE) negf/NEGFObject.h negf/selfenergies/SelfEnergy.h \
				io/Options.h negf/hamiltonian/Hamiltonian.h geometry/Overlap.h negf/GreenFunctions.h \
                negf/selfenergies/SEGolizadeh.h negf/selfenergies/SEGolizadeh.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/selfenergies/SEGolizadeh.cpp
		
obj/SEOpticalPhonon.o: $(ALL_H) $(GEOM_H) $(KSPACE_ESPACE) negf/NEGFObject.h negf/selfenergies/SelfEnergy.h \
				io/Options.h geometry/Overlap.h negf/GreenFunctions.h io/MaterialDatabase.h \
				io/PropertyContainer.h io/TernaryPropertyContainer.h \
                negf/selfenergies/SEOpticalPhonon.h negf/selfenergies/SEOpticalPhonon.cpp
	$(CPP) -c $(CPP_FLAGS) -O2 $(CPP_INC) -o $@ negf/selfenergies/SEOpticalPhonon.cpp
		
obj/SEAcousticPhonon.o: $(ALL_H) $(GEOM_H) $(KSPACE_ESPACE) negf/NEGFObject.h negf/selfenergies/SelfEnergy.h \
				io/Options.h geometry/Overlap.h negf/GreenFunctions.h io/MaterialDatabase.h \
				io/PropertyContainer.h io/TernaryPropertyContainer.h \
                negf/selfenergies/SEAcousticPhonon.h negf/selfenergies/SEAcousticPhonon.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/selfenergies/SEAcousticPhonon.cpp
		
obj/SEIonizedImpurities.o: $(ALL_H) $(GEOM_H) $(KSPACE_ESPACE) negf/NEGFObject.h \
				negf/selfenergies/SelfEnergy.h io/Options.h geometry/Overlap.h negf/GreenFunctions.h \
				io/MaterialDatabase.h io/PropertyContainer.h io/TernaryPropertyContainer.h \
                negf/selfenergies/SEIonizedImpurities.h negf/selfenergies/SEIonizedImpurities.cpp
	$(CPP) -c $(CPP_FLAGS) -O2 $(CPP_INC) -o $@ negf/selfenergies/SEIonizedImpurities.cpp
		
obj/SEPhotonSpontaneous.o: $(ALL_H) $(GEOM_H) $(KSPACE_ESPACE) negf/NEGFObject.h negf/selfenergies/SelfEnergy.h \
				io/Options.h geometry/Overlap.h negf/GreenFunctions.h io/MaterialDatabase.h \
				io/PropertyContainer.h io/TernaryPropertyContainer.h poisson/CSRMatrix.h \
				negf/hamiltonian/TdkpInfoDesk.h \
                negf/selfenergies/SEPhotonSpontaneous.h negf/selfenergies/SEPhotonSpontaneous.cpp
	$(CPP) -c $(CPP_FLAGS) -O2 $(CPP_INC) -o $@ negf/selfenergies/SEPhotonSpontaneous.cpp

obj/InnerLoop.o: $(ALL_H) $(GEOM_H) negf/*.h negf/hamiltonian/*.h negf/selfenergies/*.h \
				 negf/PostProcessing.h negf/hamiltonian/Hamiltonian.h \
                 negf/InnerLoop.h negf/InnerLoop.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/InnerLoop.cpp
	
obj/PostProcessing.o: $(ALL_H) $(GEOM_H) negf/*.h negf/hamiltonian/*.h negf/selfenergies/*.h \
				 kspace_Espace/DomainMaster.h kspace_Espace/Domain.h kspace_Espace/DomainPoint.h \
				 io/Options.h negf/selfenergies/SEPhotonSpontaneous.h negf/PostProcessing.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/PostProcessing.cpp
	
obj/Resonances.o: $(ALL_H) $(GEOM_H) negf/*.h kspace_Espace/*.h io/MaterialDatabase.h io/PropertyContainer.h \
				negf/hamiltonian/TdkpInfoDesk.h negf/Resonances.h negf/Resonances.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/Resonances.cpp
	
obj/Current.o: $(ALL_H) $(GEOM_H) negf/NEGFObject.h negf/hamiltonian/Hamiltonian.h \
				 kspace_Espace/DomainMaster.h kspace_Espace/Domain.h kspace_Espace/DomainPoint.h \
				 io/Options.h kspace_Espace/Energies.h negf/GreenFunctions.h negf/SelfEnergies.h \
                 negf/Current.h negf/Current.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/Current.cpp
	
obj/PPEmission.o: $(ALL_H) $(GEOM_H) negf/NEGFObject.h negf/hamiltonian/Hamiltonian.h \
				 kspace_Espace/DomainMaster.h kspace_Espace/Domain.h kspace_Espace/DomainPoint.h \
				 io/Options.h kspace_Espace/Energies.h negf/GreenFunctions.h negf/SelfEnergies.h \
                 negf/selfenergies/SEPhotonSpontaneous.h negf/PPEmission.h negf/PPEmission.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/PPEmission.cpp
	
obj/Luminescence.o: $(ALL_H) $(GEOM_H) negf/NEGFObject.h negf/hamiltonian/Hamiltonian.h \
				 kspace_Espace/DomainMaster.h kspace_Espace/Domain.h kspace_Espace/DomainPoint.h \
				 io/Options.h kspace_Espace/Energies.h negf/GreenFunctions.h negf/SelfEnergies.h \
                 negf/selfenergies/SEPhotonSpontaneous.h negf/Luminescence.h negf/Luminescence.cpp
	$(CPP) -c $(CPP_FLAGS) -O3 $(CPP_INC) -o $@ negf/Luminescence.cpp
	
obj/Transmission.o: $(ALL_H) $(GEOM_H) negf/NEGFObject.h \
				 kspace_Espace/DomainMaster.h kspace_Espace/Domain.h kspace_Espace/DomainPoint.h \
				 kspace_Espace/Energies.h negf/GreenFunctions.h negf/SelfEnergies.h \
                 negf/Transmission.h negf/Transmission.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/Transmission.cpp
	
obj/QuasiFermilevel.o: $(ALL_H) $(GEOM_H) negf/NEGFObject.h \
				 kspace_Espace/DomainMaster.h kspace_Espace/Domain.h kspace_Espace/DomainPoint.h \
				 kspace_Espace/Energies.h negf/GreenFunctions.h \
                 negf/QuasiFermilevel.h negf/QuasiFermilevel.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/QuasiFermilevel.cpp
		
obj/OuterLoop.o: $(ALL_H) $(GEOM_H) negf/*.h poisson/*.h \
				 negf/selfenergies/*.h negf/OuterLoop.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ negf/OuterLoop.cpp

obj/Voltages.o: $(ALL_H) $(GEOM_H) io/PropertyContainer.h io/TernaryPropertyContainer.h \
			io/InputParser.h common/Voltages.h common/Voltages.cpp
	$(CPP) -c $(CPP_FLAGS) $(CPP_INC) -o $@ common/Voltages.cpp
	
