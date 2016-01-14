/* compile using:
gcc-4.3.0 -I/home/steiger/src/release/amd64/zlib-1.2.3/include -o compression.bin compression.cc -Bdynamic -lm -L/home/steiger/src/release/amd64/lib/gcc -lstdc++ -lgcc_s -L/home/steiger/src/release/amd64/zlib-1.2.3/lib/ -lz  
*/

#include <iostream>
#include <math.h>  // for ceil
#include <cstdlib> // for exit
#include "zlib.h"

typedef unsigned long ul;

using namespace std;

int main(int argc, char **argv)
{

	int doublesize = sizeof(double);
	int charsize = sizeof(char);
	cout << "sizeof(double)=" << doublesize << ", sizeof(char)=" << charsize << endl;
	cout << "sizeof(int)=" << sizeof(int) << endl;
	
	// create a double array
	ul n = 500000;
	cout << "n=" << n << "\n\n";
	double * x = 0new double[n];
	for (ul i=0; i<n; i++) {
		x[i] = i*i;
		//cout << "x[" << i << "]=" << x[i] << ", ";
	}
	//cout << "\n\n";
	
	// cast into byte (unsigned int)
	unsigned char * xchar = (unsigned char *)(x);
	ul num_xchar = n*doublesize;
	
	// try to compress
	// from the zlib homepage: "...the destination buffer, which must be at least 0.1% larger than sourceLen plus 12 bytes."
	ul dest_size = ceil(1.01*num_xchar) + 12;
	ul * compressed_size = &dest_size;
	unsigned char * compressed = new unsigned char[dest_size];	
	cout << "uncompressed size:" << num_xchar << ", maximum compressed size:" << dest_size << "\ncompressing..."; cout.flush();
	int err = compress(compressed, compressed_size, xchar, num_xchar);
	if (err!=Z_OK) { 
		cout << "error " << err << " occurred during the compression: ";
		switch (err) {
			case Z_MEM_ERROR: cout << "out of memory."; break;
			case Z_BUF_ERROR: cout << "output buffer too small."; break;
			default: cout << "unknown.";
		}
		cout << endl; exit(1);
	}
	cout << "done. compressed size:" << *compressed_size << "\n";
	
	// imagine sending the compressed thing over the network...
	
	// at the other side:
	// from the zlib homepage: "...the destination buffer, which must be large enough to hold the entire uncompressed data."
	unsigned char * received = compressed;
	ul num_recv_chars = *compressed_size;
	ul max_size = num_xchar + 100;
	unsigned char * uncompressed = new unsigned char[max_size];
	ul * uncompressed_len = &max_size;
	
	cout << "decompressing...";cout.flush();
	err = uncompress(uncompressed, uncompressed_len, received, num_recv_chars);
	if (err!=Z_OK) { 
		cout << "error " << err << " occurred during the compression: ";
		switch (err) {
			case Z_MEM_ERROR: cout << "out of memory."; break;
			case Z_BUF_ERROR: cout << "output buffer too small."; break;
			case Z_DATA_ERROR: cout << "corrupted input data."; break;
			default: cout << "unknown.";
		}
		cout << endl; exit(1);
	}
	cout << "done. now length=" << *uncompressed_len << "\n\n";
	
	// recasting into doubles:
	double * y = (double *)(uncompressed);
	for (ul i=0; i<n; i++) {
		//cout << "y[" << i << "]=" << y[i] << ", ";
	}
	//cout << endl;
	
	// check that y and x are the same
	for (ul i=0; i<n; i++) {
		if (fabs(x[i]-y[i])>1e-14) cout << "ERROR: x["<<i<<"]="<<x[i]<<", y["<<i<<"]="<<y[i]<<"\n";
	}
	
	// clean up
	cout << "cleaning up...";
	delete [] x;
	delete [] compressed;
	delete [] uncompressed;
	cout << endl;
	
}
