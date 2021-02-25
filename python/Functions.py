
#inline double EOS(size_t const & EQ, double const & Cs0, double const & P00, double const & Density, double const & Density0)
def EOS(EQ, Cs0, P00, Density, Density0):

    if (EQ==0):
        return P00+(Cs0*Cs0)*(Density-Density0)

    else:
        if (EQ==1):
            return P00+(Density0*Cs0*Cs0/7.0)*(pow(Density/Density0,7.0)-1)
        else:
            return (Cs0*Cs0)*Density

    # default:
        # std::cout << "Please correct Pressure Equation No and run again" << std::endl
        # std::cout << "0 => P0+(Cs*Cs)*(Density-Density0)" << std::endl
        # std::cout << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)" << std::endl
        # std::cout << "2 => (Cs*Cs)*Density" << std::endl
        # abort()
        # break
    
