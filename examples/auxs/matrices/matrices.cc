#include <deal.II/lac/full_matrix.h>
#include <iostream>
#include <ostream>

int main()
{
	using namespace dealii;
	
	FullMatrix <double> B(2,4);
	FullMatrix <double> C(4,4);
	
	std::ostream out(std::cout.rdbuf() );
	
	
	B=0.0;
	B.set(0,0,1.);	B.set(0,1,2.);B.set(0,2,3.);	B.set(0,3,4.);
	B.set(1,0,5.);	B.set(1,1,6.);B.set(1,2,7.);	B.set(1,3,8.);
	
	//void 	mmult (FullMatrix< number2 > &C, const FullMatrix< number2 > &B, const bool adding=false) const
	
	// template void FullMatrix< number >::Tmmult< long double >	(	FullMatrix< number2 > & 	C,
	// const FullMatrix< number2 > & 	B,
	// const bool 	adding = false 
	// )		const

	//The optional parameter adding determines, whether the result is stored in C or added to C.
	//if (adding) C += AT*B
	//if (!adding) C = AT*B

	B.Tmmult(C,B);
	
	
	std::cout << "Out" <<std::endl;
	C.print_formatted(out);
	//out.rdbuf();
	
	// template<typename number>
	// void FullMatrix< number >::set	(	const size_type 	i,
	// const size_type 	j,
	// const number 	value 
	// )	

	// Original Array:
	// [[1 2 3 4]
	 // [5 6 7 8]]
	// Transposed Array:
	// [[1 5]
	 // [2 6]
	 // [3 7]
	 // [4 8]]
	// BtB:
	// [[26 32 38 44]
	 // [32 40 48 56]
	 // [38 48 58 68]
	 // [44 56 68 80]]

}