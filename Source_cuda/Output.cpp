#include "Domain.h"

namespace SPH {
	
inline void Domain::PrintInput(char const * FileKey)
{
	//type definition to shorten coding
	std::ostringstream oss;

	//Writing Inputs in a Log file
	String fn(FileKey);

	oss << "Dimension = "<< Dimension << "D\n";
 
	oss << "\nKernel Type = ";
	switch (KernelType)
	{
		case 0:
		oss << "Qubic Spline\n";
		break;
		case 1:
		oss << "Quintic\n";
		break;
		case 2:
		oss << "Quintic Spline\n";
		break;
	}

	oss << "\nViscosity Equation = ";
	switch (VisEq)
	{
		case 0:
			oss << "0 => Morris et al 1997\n";
			break;
		case 1:
			oss << "1 => Shao et al 2003\n";
			break;
		case 2:
			oss << "2 => Real viscosity for incompressible fluids\n";
			break;
		case 3:
			oss << "3 => Takeda et al 1994 (Real viscosity for compressible fluids)\n";
			break;
	}

	oss << "\nComputational domain size\n";
	oss << "Bottom Left-Corner Front = " << BLPF <<" m\n";
	oss << "Top Right-Corner Rear    = " << TRPR <<" m\n";

	oss << "\nMax of the smoothing lengths, h = " << hmax << " m\n";
	oss << "Cell factor in Linked List (based on kernels) = " << Cellfac << "\n";

	oss << "\nCell Size in XYZ Directions = " << CellSize <<" m\n";
	oss << "No of Cells in XYZ Directions = ( " << CellNo[0] << " , " << CellNo[1] << " , " << CellNo[2] <<" )\n" ;

	oss << "\nInitial No of Particles = " << Particles.Size() << "\n";

	oss << "\nInitial Time Step = "<<deltatint << " S\n";

	oss << "\nExternal Acceleration (Gravity enabled)= "<<Gravity<< " m/s2\n";

	oss << "\nNo of Threads = "<<Nproc<<"\n";

	oss << "\nPeriodic Boundary Condition X dir= " << (BC.Periodic[0] ? "True" : "False") << "\n";
	oss << "Periodic Boundary Condition Y dir= " << (BC.Periodic[1] ? "True" : "False") << "\n";
	oss << "Periodic Boundary Condition Z dir= " << (BC.Periodic[2] ? "True" : "False") << "\n";

	fn = FileKey;
	fn.append("_log.dat");
	std::ofstream of(fn.CStr(), std::ios::out);
	of << oss.str();
	of.close();
}

inline void Domain::WriteXDMF (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".hdf5");
    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    float * Posvec	= new float[3*Particles.Size()];
    float * Velvec	= new float[3*Particles.Size()];
    float * ACCvec	= new float[3*Particles.Size()];
    float * Pressure	= new float[  Particles.Size()];
    float * Prop1	= new float[  Particles.Size()];
    float * Prop2	= new float[  Particles.Size()];
    float * Prop3	= new float[  Particles.Size()];
    float * Density	= new float[  Particles.Size()];
    float * Mass	= new float[  Particles.Size()];
    float * sh		= new float[  Particles.Size()];
    int   * Tag			= new int  [  Particles.Size()];
    float * Sigma		= new float[6*Particles.Size()];
    float * ShearS   	= new float[6*Particles.Size()]; 	//LUCIANO
    float * Strain		= new float[6*Particles.Size()];
	float * StrainRate	= new float[6*Particles.Size()];
	float * Strain_pl	= new float[6*Particles.Size()];
	float * Sigma_eq	= new float[  Particles.Size()];	//LUCIANO
	float * Temp		= new float[  Particles.Size()];	//LUCIANO
	float * Pl_Strain	= new float[  Particles.Size()];	//LUCIANO
	int   * Nb			= new int  [  Particles.Size()];	//LUCIANO
    float * Disvec	= new float[3*Particles.Size()];		//LUCIANO
	
	double P1,P2,P3;

    #pragma omp parallel for schedule (static) private(P1,P2,P3) num_threads(Nproc)
    for (int i=0;i<Particles.Size();i++)
    {
		//LUCIANO 
		Particles[i]->CalculateEquivalentStress();
		
        Posvec  [3*i  ] = float(Particles[i]->x(0));
        Posvec  [3*i+1] = float(Particles[i]->x(1));
        Posvec  [3*i+2] = float(Particles[i]->x(2));
        Velvec  [3*i  ] = float(Particles[i]->v(0));
        Velvec  [3*i+1] = float(Particles[i]->v(1));
        Velvec  [3*i+2] = float(Particles[i]->v(2));
        ACCvec  [3*i  ] = float(Particles[i]->a(0));
        ACCvec  [3*i+1] = float(Particles[i]->a(1));
        ACCvec  [3*i+2] = float(Particles[i]->a(2));
       	Pressure[i    ] = float(Particles[i]->Pressure);
        Density [i    ] = float(Particles[i]->Density);
        Mass	[i    ] = float(Particles[i]->Mass);
        sh		[i    ] = float(Particles[i]->h);
        Tag     [i    ] = int  (Particles[i]->ID);
        Sigma   [6*i  ] = float(Particles[i]->Sigma(0,0));  //Originally0,0
        Sigma   [6*i+1] = float(Particles[i]->Sigma(1,1));	//Originally0,1
        Sigma   [6*i+2] = float(Particles[i]->Sigma(2,2));  //Originally0,2
        Sigma   [6*i+3] = float(Particles[i]->Sigma(0,1));  //Originally1,1
        Sigma   [6*i+4] = float(Particles[i]->Sigma(1,2));  //Originally1,2
        Sigma   [6*i+5] = float(Particles[i]->Sigma(0,2));  //Originally2,2
        ShearS   [6*i  ] = float(Particles[i]->ShearStress(0,0)); //Originally0,0
        ShearS   [6*i+1] = float(Particles[i]->ShearStress(1,1));//Originally0,1
        ShearS   [6*i+2] = float(Particles[i]->ShearStress(2,2)); //Originally0,2
        ShearS   [6*i+3] = float(Particles[i]->ShearStress(0,1)); //Originally1,1
        ShearS   [6*i+4] = float(Particles[i]->ShearStress(1,2)); //Originally1,2
        ShearS   [6*i+5] = float(Particles[i]->ShearStress(0,2)); //Originally2,2
        Strain  [6*i  ] = float(Particles[i]->Strain(0,0));
        Strain  [6*i+1] = float(Particles[i]->Strain(1,1));
        Strain  [6*i+2] = float(Particles[i]->Strain(2,2));
        Strain  [6*i+3] = float(Particles[i]->Strain(0,1));
        Strain  [6*i+4] = float(Particles[i]->Strain(1,2));
        Strain  [6*i+5] = float(Particles[i]->Strain(0,2));
        StrainRate  [6*i  ] = float(Particles[i]->StrainRate(0,0));
        StrainRate  [6*i+1] = float(Particles[i]->StrainRate(1,1));
        StrainRate  [6*i+2] = float(Particles[i]->StrainRate(2,2));
        StrainRate  [6*i+3] = float(Particles[i]->StrainRate(0,1));
        StrainRate  [6*i+4] = float(Particles[i]->StrainRate(1,2));
        StrainRate  [6*i+5] = float(Particles[i]->StrainRate(0,2));
        Strain_pl  [6*i  ] = float(Particles[i]->Strain_pl(0,0));
        Strain_pl  [6*i+1] = float(Particles[i]->Strain_pl(1,1));
        Strain_pl  [6*i+2] = float(Particles[i]->Strain_pl(2,2));
        Strain_pl  [6*i+3] = float(Particles[i]->Strain_pl(0,1));
        Strain_pl  [6*i+4] = float(Particles[i]->Strain_pl(1,2));
        Strain_pl  [6*i+5] = float(Particles[i]->Strain_pl(0,2));
		Sigma_eq[i    ] = float(Particles[i]->Sigma_eq);
		Temp	[i    ] = float(Particles[i]->T);
		Pl_Strain	[i] = float(Particles[i]->pl_strain);
		Nb		[i    ] = int(Particles[i]->Nb);
		
		
        Disvec  [3*i  ] = float(Particles[i]->Displacement(0));
        Disvec  [3*i+1] = float(Particles[i]->Displacement(1));
        Disvec  [3*i+2] = float(Particles[i]->Displacement(2));		
		
	UserOutput(Particles[i],P1,P2,P3);
        Prop1	[i    ] = float(P1);
        Prop2	[i    ] = float(P2);
        Prop3	[i    ] = float(P3);
   }

    int data[1];
    String dsname;
    hsize_t dims[1];
    dims[0]=1;
    data[0]=Particles.Size();
    dsname.Printf("/NP");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,data);
    dims[0] = 3*Particles.Size();
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("Velocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dsname.Printf("Acceleration");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,ACCvec);
    dims[0] = Particles.Size();
    dsname.Printf("Tag");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tag);
    dsname.Printf("Pressure");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Pressure);
    dsname.Printf("Density");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Density);
    dsname.Printf(OutputName[0]);
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Prop1);
    dsname.Printf(OutputName[1]);
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Prop2);
    dsname.Printf(OutputName[2]);
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Prop3);
    dsname.Printf("Mass");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Mass);
    dsname.Printf("h");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,sh);
    dims[0] = 6*Particles.Size();
    dsname.Printf("Sigma");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Sigma);
    dsname.Printf("ShearS");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Sigma);
    dsname.Printf("Strain");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Strain);
    dsname.Printf("StrainRate");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,StrainRate);
    dsname.Printf("Strain_pl");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Strain_pl);
    dims[0] = Particles.Size();
	dsname.Printf("Sigma_eq");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Sigma_eq);
	dsname.Printf("Temperature");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Temp);
    dsname.Printf("Pl_Strain");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Pl_Strain);
    dsname.Printf("Neighbors");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Nb);
    dims[0] = 3*Particles.Size();
	dsname.Printf("Displacement");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Disvec);
	
    delete [] Posvec;
    delete [] Velvec;
    delete [] ACCvec;
    delete [] Pressure;
    delete [] Prop1;
    delete [] Prop2;
    delete [] Prop3;
    delete [] Density;
    delete [] Mass;
    delete [] sh;
    delete [] Tag;
    delete [] Sigma;
    delete [] ShearS;
    delete [] Strain;
    delete [] StrainRate;
    delete [] Strain_pl;
	delete [] Temp;
	delete [] Sigma_eq;
	delete [] Pl_Strain;
	delete [] Nb;
	delete [] Disvec;
	
   //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"SPHCenter\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"10\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Position\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Acceleration\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Acceleration \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Density\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Density \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Pressure \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"" << OutputName[0] << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/" << OutputName[0] << " \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"" << OutputName[1] << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/" << OutputName[1] << " \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"" << OutputName[2] << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/" << OutputName[2] << " \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Sigma\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Sigma \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
	oss << "     <Attribute Name=\"ShearS\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/ShearS \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Strain\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Strain \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"StrainRate\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/StrainRate \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Strain_pl\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Strain_pl \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Temperature\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Temperature \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Sigma_eq\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Sigma_eq \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Pl_Strain\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Pl_Strain \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Neighbors\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Neighbors \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Displacement\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Displacement \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";


    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

}; // namespace SPH