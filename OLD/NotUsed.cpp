inline void Domain::AddBoxNo(int tag, Vec3_t const & V, size_t nx, size_t ny, size_t nz, double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
{
    if (!(type==0 || type==1))
    {
	   	std::cout << "Packing Type is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Hexagonal Close Packing" << std::endl;
		std::cout << "1 => Cubic Packing" << std::endl;
	    abort();
    }

    if (!(rotation==0 || rotation==90))
    {
	   	std::cout << "Packing Rotation Angle is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => " << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << std::endl;
		std::cout << "90 =>" << std::endl;
		std::cout << "  0   0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0   0  " << std::endl;
		abort();
    }

//	Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating particles by AddBoxNo with defined numbers of particles--------------" << std::endl;

    size_t PrePS = Particles.Size();

    double x,y;

    double qin = 0.03;
    srand(100);

    if (Dimension==3)
    {
   		double z;

    	if (type==0)
    	{
    		//Hexagonal close packing
 		    for (size_t k=0; k<nz; k++)
		    for (size_t j=0; j<ny; j++)
		    for (size_t i=0; i<nx; i++)
		    {
				if ((k%2!=0) && (j%2!=0)) x = V(0) + (2*i+(j%2)+(k%2)-1)*r; else x = V(0) + (2*i+(j%2)+(k%2)+1)*r;
				y = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
				z = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
				if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					else    Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
			}
    	}
    	else
    	{
    		//Cubic packing
		    for (size_t k=0; k<nz; k++)
		    for (size_t j=0; j<ny; j++)
		    for (size_t i=0; i<nx; i++)
		    {
				x = V(0) + (2.0*i+1)*r;
				y = V(1) + (2.0*j+1)*r;
				z = V(2) + (2.0*k+1)*r;
				if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					else    Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
			}
    	}

        //Calculate particles' mass in 3D
        Vec3_t temp, Max=V;
		for (size_t i=PrePS; i<Particles.Size(); i++)
		{
			if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		}
		Max +=r;
		temp = Max-V;
		double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.Size()-PrePS);

		#pragma omp parallel for num_threads(Nproc)
		for (size_t i=PrePS; i<Particles.Size(); i++)
		{
			Particles[i]->Mass = Mass;
		}
    }

    if (Dimension==2)
    {
    	if (type==0)
    	{
    		//Hexagonal close packing
    		if (rotation==0)
    		{
    		    for (size_t j=0; j<ny; j++)
    		    for (size_t i=0; i<nx; i++)
    		    {
					x = V(0) + (2*i+(j%2)+1)*r;
					y = V(1) + (sqrt(3.0)*j+1)*r;
					if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
				}
			}
    		else
    		{
    		    for (size_t i=0; i<nx; i++)
    		    for (size_t j=0; j<ny; j++)
    		    {
					x = V(0) + (sqrt(3.0)*i+1)*r;
					y = V(1) + (2*j+(i%2)+1)*r;
					if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
				}
    		}
    	}
    	else
    	{
    		//Cubic packing
		    for (size_t j=0; j<ny; j++)
		    for (size_t i=0; i<nx; i++)
		    {
				x = V(0) + (2*i+1)*r;
				y = V(1) + (2*j+1)*r;
				if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
					else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),2.0*r*2.0*r*Density,Density,h,Fixed));
			}

    	}
    }

	R = r;
}