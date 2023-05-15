#include "Domain.h"

namespace SPH {
void Domain::ReadXDMF			(char const * FileKey) {
    
    hid_t file_id;
    hid_t acc_tpl1;		/* File access templates */    
    file_id=H5Fopen(FileKey,H5F_ACC_RDWR,acc_tpl1);
    assert(fid1 != FAIL);
    //MESG("H5Fopen succeed");
    
    //file_id = H5Fopen(name, flags = h5default("H5F_ACC_RD"), fapl = NULL, native = FALSE)
    
    int data[1];
    String dsname;
    hsize_t dims[1];
    dims[0]=1;
    //data[0]=Particles.Size();
    dsname.Printf("/NP");

    H5LTread_dataset_int(file_id,dsname.CStr(),data);
    
    cout << "Particle Size "<<data[0]<<endl;
    
    dims[0] = 3*Particles.Size();
    dsname.Printf("Position");
    
    float * Posvec	= new float[3*data[0]];
    H5LTread_dataset_float ( file_id, dsname.CStr(),Posvec);
    for (int i=0;i<data[0];i++){
      //cout << Posvec[3*i]<<", "<<Posvec[3*i+1]<<", "<<Posvec[3*i+2]<<endl;
      //Particles.Push(new Particle(tag,Vec3_t(Posvec[3*i],Posvec[3*i+1],Posvec[3*i+2]),Vec3_t(0,0,0),0.0,Density,h,false));	
      Particles.Push(new Particle(0,Vec3_t(Posvec[3*i],Posvec[3*i+1],Posvec[3*i+2]),Vec3_t(0,0,0),0.0,1.0,1.0,false));		
    }
    // H5LTread_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    // cout << "Buffer Size"<<endl;
    
    // dsname.Printf("Velocity");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    // dsname.Printf("Acceleration");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,ACCvec);
    // dims[0] = Particles.Size();
    // dsname.Printf("Tag");
    // H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tag);
    // dsname.Printf("Pressure");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Pressure);
    // dsname.Printf("Density");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Density);
    // dsname.Printf(OutputName[0]);
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Prop1);
    // dsname.Printf(OutputName[1]);
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Prop2);
    // dsname.Printf(OutputName[2]);
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Prop3);
    // dsname.Printf("Mass");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Mass);
    // dsname.Printf("h");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,sh);
    // dims[0] = 6*Particles.Size();
    // dsname.Printf("Sigma");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Sigma);
    // dsname.Printf("ShearS");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Sigma);
    // dsname.Printf("Strain");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Strain);
    // dsname.Printf("StrainRate");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,StrainRate);
    // dsname.Printf("Strain_pl");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Strain_pl);
    // dsname.Printf("gradcorrmat");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,gradcorrmat);
    // dims[0] = Particles.Size();
	// dsname.Printf("Sigma_eq");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Sigma_eq);
	// dsname.Printf("Temperature");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Temp);
    // dsname.Printf("Pl_Strain");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Pl_Strain);
    // dsname.Printf("Neighbors");
    // H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Nb);
    // dsname.Printf("ContNeib");
    // H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,ContNb);
    // dsname.Printf("q_friction");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,q_friction);
    // dsname.Printf("c_shearabs");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,c_shearabs);

    // dsname.Printf("ps_en");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,ps_en);
    
    // dims[0] = 3*Particles.Size();
	// dsname.Printf("Displacement");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Disvec);
	// dsname.Printf("Contact Force");
    // H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,ContForce);	
  
}

};