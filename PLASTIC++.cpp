//This is the main code to run PLASTIC in C++.  Go to line 121 to input the simulation parameters.

#include <iostream>
#include <cmath>
#include <vector>
#include <fftw3.h>
#include <chrono>
#include <complex>
#include <fstream>
#include <string>
#include <algorithm>


using namespace std;
const double pi = 3.14159265358979323846;



int indexWrap(int index, int size) {
    // Wrap negative indices to the valid range
    return (index % size + size) % size;
}




std::string removeChars(const std::string& str, const std::string& charsToRemove) {
    std::string result;
    for (char c : str) {
        if (charsToRemove.find(c) == std::string::npos) {
            result += c;
        }
    }
    return result;
}

// Function to split a string into a vector of strings based on a delimiter
std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::vector<std::vector<std::complex<double>>>> create3DVector(int meshpoints) {
    return std::vector<std::vector<std::vector<std::complex<double>>>>(meshpoints, 
        std::vector<std::vector<std::complex<double>>>(meshpoints, 
        std::vector<std::complex<double>>(meshpoints, std::complex<double>(0.0, 0.0))));
}
//vector<vector<vector<complex<double>>>> create3DVector(int size) 
//{
//    return vector<vector<vector<complex<double>>>>(size, vector<vector<complex<double>>>(size, vector<complex<double>>(size)));
//}
vector<complex<double>> create1DVector(int size)
{
    return vector<complex<double>>(size);
}
vector<double> create1DVectorReal(int size)
{
    return vector<double>(size);
}


// Flatten function for 3D vector of doubles to 1D vector of complex doubles
std::vector<std::complex<double>> flatten(int meshpoints, const std::vector<std::vector<std::vector<double>>>& integrand) {
    std::vector<std::complex<double>> integrand_flat(meshpoints * meshpoints * meshpoints);
    for (int i = 0; i < meshpoints; i++) {
        for (int j = 0; j < meshpoints; j++) {
            for (int k = 0; k < meshpoints; k++) {
                int index = i * meshpoints * meshpoints + j * meshpoints + k;
                integrand_flat[index] = std::complex<double>(integrand[i][j][k], 0.0);
            }
        }
    }
    return integrand_flat;
}

// Flatten function for 3D vector of complex doubles to 1D vector of complex doubles
std::vector<std::complex<double>> flatten(int meshpoints, const std::vector<std::vector<std::vector<std::complex<double>>>>& integrand) {
    std::vector<std::complex<double>> integrand_flat(meshpoints * meshpoints * meshpoints);
    for (int i = 0; i < meshpoints; i++) {
        for (int j = 0; j < meshpoints; j++) {
            for (int k = 0; k < meshpoints; k++) {
	      int index = (i * meshpoints * meshpoints + j * meshpoints + k);
                integrand_flat[index] = integrand[i][j][k];
            }
        }
    }
    return integrand_flat;
}

// Unflatten function to convert 1D vector of complex doubles back to 3D vector of complex doubles
std::vector<std::vector<std::vector<std::complex<double>>>> unflatten(int meshpoints, const std::vector<std::complex<double>>& integrand_flat) {
    std::vector<std::vector<std::vector<std::complex<double>>>> integrand(meshpoints, std::vector<std::vector<std::complex<double>>>(meshpoints, std::vector<std::complex<double>>(meshpoints)));
    for (int i = 0; i < meshpoints; ++i) {
        for (int j = 0; j < meshpoints; ++j) {
            for (int k = 0; k < meshpoints; ++k) {
                int index = i * meshpoints * meshpoints + j * meshpoints + k;
                integrand[i][j][k] = integrand_flat[index];
            }
        }
    }
    return integrand;
}void print_complex_array(fftw_complex* arr, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << "(" << arr[i][0] << ", " << arr[i][1] << ") ";
    }
    std::cout << std::endl;

}

int main()
{

    double Omega_m;
    double H0;
    double z_init;
    double a_init;
    double z_final;
    double a_final;
    double boxsize;
    int meshpoints;
    double sigma;
    double cs2_0;
    double cv2_0;
    string initial_conditions_path;
    int auto_stop_enabled;
    int max_steps;
    double step_factor;
    string pion_output_path;
    string pion_prime_output_path;
    string overdensity_output_path;
    

    string line;
    string var;
    string val;
    ifstream inputsfile("inputs.txt");
    if (inputsfile.is_open()) {
        while (getline(inputsfile, line)) {
	  if (line[0] != '#' and line.length() > 1){
	    line = line.substr(0, line.find("#"));
	    var = removeChars(line.substr(0, line.find("="))," ");
	    val = removeChars(line.substr(line.find("=") + 1,line.length()), " ");
	    if (var == "Omega_m") {
	      Omega_m = std::stod(val);
	    }
	    if (var == "H0") {
	      H0 = std::stod(val);
	    }
	    if (var == "z_init") {
	      z_init = std::stod(val);
	      a_init = 1.0/(1.0 + z_init);
	    }
	    if (var == "a_init") {
	      a_init = std::stod(val);
	      z_init = 1.0/a_init - 1.0;
	    }
	    if (var == "z_final") {
	      z_final = std::stod(val);
	      a_final = 1.0/(1.0 + z_final);
	    }
	    if (var == "a_final") {
	      a_final = std::stod(val);
	      z_final = 1.0/a_final - 1.0;
	    }
	    if (var == "boxsize") {
	      boxsize = std::stod(val);
	    }
	    if (var == "meshpoints") {
	      meshpoints = std::stoi(val);
	    }
	    if (var == "sigma") {
	      sigma = std::stod(val);
	    }
	    if (var == "cs2_0") {
	      cs2_0 = std::stod(val);
	    }
	    if (var == "cv2_0") {
	      cv2_0 = std::stod(val);
	    }
	    if (var == "initial_conditions_path") {
	      initial_conditions_path = val;
	    }
	    if (var == "auto_stop_enabled") {
	      auto_stop_enabled = std::stoi(val);
	    }
	    if (var == "max_steps") {
	      max_steps = std::stoi(val);
	    }
	    if (var == "step_factor") {
	      step_factor = std::stod(val);
	    }
	    if (var == "pion_output_path") {
	      pion_output_path = val;
	    }
	    if (var == "pion_prime_output_path") {
	      pion_prime_output_path = val;
	    }
	    if (var == "overdensity_output_path") {
	      overdensity_output_path = val;
	    }
	  }
        }
        inputsfile.close();
    }
  // //Input cosmic and simulation parameters.  The first section should be identical to the parameters
  // //in the initial conditions file.

  //   //Cosmic
  //   double Omega_m = 0.31; // matter density fraction in present epoch
  //   double H0 = 68.0; // km/s/Mpc
  //   double z_init = 99.0; // initial redshift
    
  //   //Simulation
  //   double boxsize = 600.0; // Mpc
  //   const int meshpoints = 20;

  //   //Input and output file names.
  //   std::string initial_conditions_path = "initial_conditions.txt";
  //   std::string pion_output_path = "PidataC++.txt";
  //   std::string pion_prime_output_path = "PiPrimedataC++.txt";
  //   std::string overdensity_output_path = "DeltaC++.txt";
    
  //   double sigma = 25.0; // Smoothing scale for Gaussian filter in Mpc.
  //   double z_final = 0;  // Final redshift
  //   double cs2_0 = 0E-6; // Time dependence assumed to be D^2
  //   double cv2_0 = 0E-6;  // Timde dependence assumed to be D^2/sqrt(Omega_M)
  //   int max_steps = 6000;


    
    int total_points = meshpoints * meshpoints * meshpoints;
    double meshsize = boxsize / meshpoints;
    std::cout << "Size of periodic box is: " << boxsize << " Mpc" << std::endl;
    std::cout << "Number of grid points is: " << meshpoints << std::endl;
    cout << "Initial conditions at z = " << z_init << " read from " << initial_conditions_path << "\n";
    

    a_init = 1.0 / (1.0 + z_init);
    a_final = 1.0 / (1.0 + z_final);
    double Hubble_init = (H0 / 299792.458) * a_init * std::sqrt(Omega_m / std::pow(a_init, 3));
    double time_init = 2.0 * std::sqrt(a_init) * 299792.458 / (H0 * std::sqrt(Omega_m));
    double time = time_init;
        
    double mesh[meshpoints][meshpoints][meshpoints][3] = { 0.0 };
    cout << "Type OK to continue:";
    int abc; cin >> abc;
    
    string filePath = initial_conditions_path;
    
    std::vector<std::vector<std::vector<double>>> ic1(
        meshpoints, std::vector<std::vector<double>>(
            meshpoints, std::vector<double>(meshpoints, 0.0)
        )
    );

    std::vector<std::vector<std::vector<double>>> icprime1(
        meshpoints, std::vector<std::vector<double>>(
            meshpoints, std::vector<double>(meshpoints, 0.0)
        )
    );

    auto ic = flatten(meshpoints, ic1);
    auto icprime = flatten(meshpoints, icprime1);

    std::ifstream file(filePath);
    if (!file) {
        std::cerr << "Unable to open file " << filePath << std::endl;
        return 1;
    }
    else {
        std::cout << "file opened" << std::endl;
    }

    std::string fileContent((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    std::string flattened = removeChars(removeChars(fileContent,"["),"]");

    std::vector<std::string> ic_strings = split(flattened, ',') ;

    for (int i = 0; i < meshpoints * meshpoints * meshpoints; ++i){
      ic[i] = std::stod(ic_strings[i]);
      icprime[i] = ic[i]/time_init;
    }
    
    auto pidata = create1DVector(meshpoints* meshpoints* meshpoints);
    auto pidatanext = create1DVector(total_points);

    vector<complex<double>> piprimedata(total_points, 0.0);
    vector<complex<double>> piprimedatanext(total_points, 0.0);
    pidata = ic;
    pidatanext = pidata;
    piprimedata = icprime;
    piprimedatanext = icprime;
    
    time = time_init;
    double a_scale = a_init;
    double lin_pi_growth_times_a2 = ( (2.0 * sqrt( pow(a_init, 5) )* 299792.458) / (5.0 * H0 *sqrt(Omega_m)));

    for (int i = 1; i < 1000000; i++){
      if (a_scale > 1.0){
	break;
	}
      double a_step = a_scale * 0.0001;
      lin_pi_growth_times_a2 += ((a_step * 299792.458) / (H0 * sqrt((Omega_m/pow(a_scale,3))+(1-Omega_m)) ) );
      a_scale += a_step;
    }

    double D_zero = lin_pi_growth_times_a2;

    a_scale = a_init;
    lin_pi_growth_times_a2 = ((2.0 * sqrt( pow(a_init, 5) ) *  299792.458) / (5.0 * H0 *sqrt(Omega_m))  );
    
    
    int numsteps = 0;

    auto overdensity = create1DVector(total_points);
    auto theta = create1DVector(total_points);
    auto thetaprime = create1DVector(total_points);
    auto delta = create1DVector(total_points);
    vector<complex<double>> NewtonPhi(total_points, 0.0);
    vector<complex<double>> NewtonPhi_FT(total_points, 0.0);
    vector<complex<double>> pi_FT(total_points, 0.0);
    vector<complex<double>> delta_FT(total_points, 0.0);
    vector<complex<double>> theta_FT(total_points, 0.0);
    vector<complex<double>> thetaprime_FT(total_points, 0.0);
    vector<complex<double>> pi_prime_FT(total_points, 0.0);
    vector<complex<double>> vx_FT(total_points, 0.0);
    vector<complex<double>> vy_FT(total_points, 0.0);
    vector<complex<double>> vz_FT(total_points, 0.0);
    vector<complex<double>> vx(total_points, 0.0);
    vector<complex<double>> vy(total_points, 0.0);
    vector<complex<double>> vz(total_points, 0.0);
    vector<complex<double>> vx_prime_FT(total_points, 0.0);
    vector<complex<double>> vy_prime_FT(total_points, 0.0);
    vector<complex<double>> vz_prime_FT(total_points, 0.0);
    vector<complex<double>> vx_prime(total_points, 0.0);
    vector<complex<double>> vy_prime(total_points, 0.0);
    vector<complex<double>> vz_prime(total_points, 0.0);
    vector<complex<double>> integrand(total_points, 0.0);
    vector<complex<double>> FTinverseLaplacian(total_points, 0.0);
    vector<complex<double>> integral(total_points, 0.0);
    vector<complex<double>> gradsquared(total_points, 0.0);
    vector<complex<double>> gradpipiprime(total_points, 0.0);
    vector<complex<double>> gradx_delta_FT(total_points, 0.0);
    vector<complex<double>> grady_delta_FT(total_points, 0.0);
    vector<complex<double>> gradz_delta_FT(total_points, 0.0);
    vector<complex<double>> gradx_delta(total_points, 0.0);
    vector<complex<double>> grady_delta(total_points, 0.0);
    vector<complex<double>> gradz_delta(total_points, 0.0);
    vector<complex<double>> piprimedatanext_FT(total_points, 0.0);
    vector<complex<double>> resetting(total_points, 0.0);
    

    auto start = chrono::high_resolution_clock::now();
    // std::cout << "The code still works at this point. " << std::endl;
    for (int t = 0; t < max_steps; t++)
    {
        
        if (a_scale > a_final) {
            break;
        }
        double a_step = a_scale * step_factor;
        double timestep = (a_step * 299792.458) / ((H0) * pow(a_scale,2) * sqrt((Omega_m / pow(a_scale,3)) + (1 - Omega_m)));
        double Hubble = (H0 / 299792.458) * a_scale * sqrt((Omega_m / pow(a_scale, 3)) + (1 - Omega_m));
	double Omega_M = (Omega_m / (Omega_m + (1.0 - Omega_m)*pow(a_scale,3)));
        double Hubble_prime = (1.0 - ((3.0 / 2.0) * Omega_M)) * pow(Hubble,2);
	double D_growth = (lin_pi_growth_times_a2 / pow(a_scale,2));
	double tildec2_s = (1.5 * Omega_M - 3.5)*(pow(D_growth,2)/(pow(D_zero,2)))*sqrt(Omega_m/Omega_M)*cv2_0;
       	tildec2_s += ((2.0 * D_growth)/(Hubble * pow(D_zero,2)))*sqrt(Omega_m/Omega_M)*cv2_0 + (pow(D_growth,2)/pow(D_zero,2))*cs2_0;
	tildec2_s += (2.0/(3.0*Omega_M))*(-3.0*(pow(D_growth,2)/pow(D_zero,2))*cs2_0 + ((2.0 * D_growth)/(Hubble * pow(D_zero,2)))*cs2_0);
	double tildec2_v_overH = (1.0/Hubble)*pow(D_growth/D_zero,2)*sqrt(Omega_m/Omega_M)*cv2_0;
	tildec2_v_overH += (1/Hubble)*(2.0/(3.0*Omega_M))*(-3.0*pow(D_growth/D_zero,2)+ ((2.0 * D_growth)/(Hubble * pow(D_zero,2))))*cs2_0;
	
	NewtonPhi = resetting;
	gradsquared = resetting;
	integrand = resetting;
	integral = resetting;
	piprimedatanext = resetting;
	piprimedatanext_FT = resetting;

        fftw_plan plan_forward_pi = fftw_plan_dft_3d(
            meshpoints,
            meshpoints,
            meshpoints,
            reinterpret_cast<fftw_complex*>(pidata.data()),
            reinterpret_cast<fftw_complex*>(pi_FT.data()),
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );
	
        fftw_plan plan_forward_pi_prime = fftw_plan_dft_3d(
            meshpoints,
            meshpoints,
            meshpoints,
            reinterpret_cast<fftw_complex*>(piprimedata.data()),
            reinterpret_cast<fftw_complex*>(pi_prime_FT.data()),
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );

        fftw_execute(plan_forward_pi);
        fftw_execute(plan_forward_pi_prime);

	//Gaussian smoothing

	double k_cutoff = 6.28318/(sigma + 0.00001);

        for (int i = 0; i < meshpoints; ++i) {
            for (int j = 0; j < meshpoints; ++j) {
                for (int k = 0; k < meshpoints; ++k) {
                    double k_mag = (6.28318 / boxsize) * std::sqrt(std::pow(std::min(i, meshpoints - i), 2) +
                        std::pow(std::min(j, meshpoints - j), 2) +
                        std::pow(std::min(k, meshpoints - k), 2));
                    int idx = i * meshpoints * meshpoints + j * meshpoints + k;
                    pi_FT[idx] *= exp( - pow(k_mag, 2)/pow(k_cutoff, 2) );
                    pi_prime_FT[idx] *= exp( - pow(k_mag, 2)/pow(k_cutoff, 2) );
                }
            }
        }
        fftw_plan plan_inverse_pi = fftw_plan_dft_3d(
            meshpoints,
            meshpoints,
            meshpoints,
            reinterpret_cast<fftw_complex*>(pi_FT.data()),
            reinterpret_cast<fftw_complex*>(pidata.data()),
            FFTW_BACKWARD,
            FFTW_ESTIMATE
        );
        fftw_execute(plan_inverse_pi);

        // Normalize and extract the real part of the result into pidata
        double norm_factor = 1.0 / (meshpoints * meshpoints * meshpoints );
        for (int i = 0; i < meshpoints; ++i) {
            for (int j = 0; j < meshpoints; ++j) {
                for (int k = 0; k < meshpoints; ++k) {
                    int idx = i * meshpoints * meshpoints + j * meshpoints + k;
                    pidata[idx] = pidata[idx].real() * norm_factor;
                }
            }
        }

        fftw_plan plan_inverse_pi_prime = fftw_plan_dft_3d(
            meshpoints,
            meshpoints,
            meshpoints,
            reinterpret_cast<fftw_complex*>(pi_prime_FT.data()),
            reinterpret_cast<fftw_complex*>(piprimedata.data()),
            FFTW_BACKWARD,
            FFTW_ESTIMATE
        );

        // Execute the inverse FFT plan
        fftw_execute(plan_inverse_pi_prime);
        // Normalize and extract the real part of the result into pidatanext
        for (int i = 0; i < meshpoints; ++i) {
            for (int j = 0; j < meshpoints; ++j) {
                for (int k = 0; k < meshpoints; ++k) {
                    int idx = i * meshpoints * meshpoints + j * meshpoints + k;
                    piprimedata[idx] = piprimedata[idx].real() * norm_factor ;
                }
            }
        }

        
        // Clean up the FFTW plan
        fftw_destroy_plan(plan_inverse_pi);
        fftw_destroy_plan(plan_inverse_pi_prime);
        fftw_destroy_plan(plan_forward_pi);
        fftw_destroy_plan(plan_forward_pi_prime);
	fftw_cleanup();


        //int term = 0;
        for (int i = 0; i < meshpoints; i++) {
	  double kx = (6.28318/boxsize)*double(i);
	  if (i > meshpoints/2) {
	    kx -= (6.28318/boxsize)*double(meshpoints);
	  }
            for (int j = 0; j < meshpoints; j++) {
	      	  double ky = (6.28318/boxsize)*double(j);
		  if (j > meshpoints/2) {
		    ky -= (6.28318/boxsize)*double(meshpoints);
		  }
                for (int k = 0; k < meshpoints; k++) {
		  double kz = (6.28318/boxsize)*double(k);
		  if (k > meshpoints/2) {
		    kz -= (6.28318/boxsize)*double(meshpoints);
		  }
                    // Flattened index for current point
                    int idx = i * meshpoints * meshpoints + j * meshpoints + k;
		    double k_mag = std::sqrt(kx * kx + ky * ky + kz * kz);
		    complex <double> imag_unit (0.0, 1.0);

		    vx_FT[idx] = - imag_unit * kx * pi_FT[idx];
		    vy_FT[idx] = - imag_unit * ky * pi_FT[idx];
		    vz_FT[idx] = - imag_unit * kz * pi_FT[idx];
		    vx_prime_FT[idx] = - imag_unit * kx * pi_prime_FT[idx];
		    vy_prime_FT[idx] = - imag_unit * ky * pi_prime_FT[idx];
		    vz_prime_FT[idx] = - imag_unit * kz * pi_prime_FT[idx];

                }
            }
        }

	fftw_plan plan_backward_vx = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(vx_FT.data()), reinterpret_cast<fftw_complex*>(vx.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_vx);

	fftw_plan plan_backward_vy = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(vy_FT.data()), reinterpret_cast<fftw_complex*>(vy.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_vy);

	fftw_plan plan_backward_vz = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(vz_FT.data()), reinterpret_cast<fftw_complex*>(vz.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_vz);
	
	fftw_plan plan_backward_vx_prime = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(vx_prime_FT.data()), reinterpret_cast<fftw_complex*>(vx_prime.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_vx_prime);

	fftw_plan plan_backward_vy_prime = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(vy_prime_FT.data()), reinterpret_cast<fftw_complex*>(vy_prime.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_vy_prime);

	fftw_plan plan_backward_vz_prime = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(vz_prime_FT.data()), reinterpret_cast<fftw_complex*>(vz_prime.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_vz_prime);	
       

	fftw_destroy_plan(plan_backward_vx);
	fftw_destroy_plan(plan_backward_vy);
	fftw_destroy_plan(plan_backward_vz);
	fftw_destroy_plan(plan_backward_vx_prime);
	fftw_destroy_plan(plan_backward_vy_prime);
	fftw_destroy_plan(plan_backward_vz_prime);
        fftw_cleanup();

	for (int i = 0; i < meshpoints; ++i) {
            for (int j = 0; j < meshpoints; ++j) {
                for (int k = 0; k < meshpoints; ++k) {
                    int idx = i * meshpoints * meshpoints + j * meshpoints + k;
                    vx[idx] = vx[idx].real() * norm_factor;
		    vy[idx] = vy[idx].real() * norm_factor;
		    vz[idx] = vz[idx].real() * norm_factor;
		    vx_prime[idx] = vx_prime[idx].real() * norm_factor;
		    vy_prime[idx] = vy_prime[idx].real() * norm_factor;
		    vz_prime[idx] = vz_prime[idx].real() * norm_factor;
		    gradsquared[idx] = vx[idx] * vx[idx] + vy[idx] * vy[idx] + vz[idx] * vz[idx];
		    gradpipiprime[idx] = vx[idx] * vx_prime[idx] + vy[idx] * vy_prime[idx] + vz[idx] * vz_prime[idx];

		    NewtonPhi[idx] = -(piprimedata[idx] + (Hubble)*pidata[idx] + (1.0 / 2.0) * gradsquared[idx]);
		}
	    }
	}

	fftw_plan plan_forward_NewtonPhi = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(NewtonPhi.data()), reinterpret_cast<fftw_complex*>(NewtonPhi_FT.data()), FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan_forward_NewtonPhi);

	for (int i = 0; i < meshpoints; ++i) {
            for (int j = 0; j < meshpoints; ++j) {
                for (int k = 0; k < meshpoints; ++k) {
		    int idx = i * meshpoints * meshpoints + j * meshpoints + k;
		    double factor = 6.28318 * 6.28318;
		    int idx_i = std::min(i, meshpoints - i);
                    int idx_j = std::min(j, meshpoints - j);
                    int idx_k = std::min(k, meshpoints - k);
		    double k_squared = double(idx_i * idx_i + idx_j * idx_j + idx_k * idx_k);

                    delta_FT[idx] = (NewtonPhi_FT[idx] * factor * k_squared) / ( boxsize * boxsize);
                    theta_FT[idx] = (-pi_FT[idx] * factor * k_squared) / (boxsize * boxsize);
		    thetaprime_FT[idx] = (-pi_prime_FT[idx] * factor * k_squared) / (boxsize * boxsize);
                }
            }
        }
        
        fftw_plan plan_backward_delta = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(delta_FT.data()), reinterpret_cast<fftw_complex*>(delta.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_delta);
              
        fftw_plan plan_backward_theta = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(theta_FT.data()), reinterpret_cast<fftw_complex*>(theta.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_theta);

	fftw_plan plan_backward_thetaprime = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(thetaprime_FT.data()), reinterpret_cast<fftw_complex*>(thetaprime.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_thetaprime);

        fftw_destroy_plan(plan_forward_NewtonPhi);
        fftw_destroy_plan(plan_backward_delta);
        fftw_destroy_plan(plan_backward_theta);
	fftw_destroy_plan(plan_backward_thetaprime);
        fftw_cleanup();

	bool autostop = false;
	
        for (int i = 0; i < meshpoints; i++) {
            for (int j = 0; j < meshpoints; j++) {
                for (int k = 0; k < meshpoints; k++) {
                    int idx = i * meshpoints * meshpoints + j * meshpoints + k;
                    overdensity[idx] = delta[idx].real() * norm_factor / (Hubble_prime - Hubble * Hubble);
                    if (overdensity[idx].real() <-1.0){ 
			autostop = true;
                    }
                    delta[idx] = delta[idx].real() * norm_factor;
                    theta[idx] *= norm_factor;
		    thetaprime[idx] *= norm_factor;
                }
            }
        }

	if (autostop == true){
	  cout << "Underdensity is less than -1.0 at time " << time << " Mpc\n";
	  //cout << "Max overdensity is " << std::max_element(overdensity);
	  cout << "Stopping the run.";
	  break;
	}

        for (int i = 0; i < meshpoints; i++) {
	  double kx = (6.28318/boxsize)*double(i);
	  if (i > meshpoints/2) {
	    kx -= (6.28318/boxsize)*double(meshpoints);
	  }
            for (int j = 0; j < meshpoints; j++) {
	      	  double ky = (6.28318/boxsize)*double(j);
		  if (j > meshpoints/2) {
		    ky -= (6.28318/boxsize)*double(meshpoints);
		  }
                for (int k = 0; k < meshpoints; k++) {
		  double kz = (6.28318/boxsize)*double(k);
		  if (k > meshpoints/2) {
		    kz -= (6.28318/boxsize)*double(meshpoints);
		  }
                    // Flattened index for current point
                    int idx = i * meshpoints * meshpoints + j * meshpoints + k;
		    double k_mag = std::sqrt(kx * kx + ky * ky + kz * kz);
		    complex <double> imag_unit (0.0, 1.0);

		    gradx_delta_FT[idx] = - imag_unit * kx * delta_FT[idx];
		    grady_delta_FT[idx] = - imag_unit * ky * delta_FT[idx];
		    gradz_delta_FT[idx] = - imag_unit * kz * delta_FT[idx];

                }
            }
        }

	fftw_plan plan_backward_gradx_delta = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(gradx_delta_FT.data()), reinterpret_cast<fftw_complex*>(gradx_delta.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_gradx_delta);
              
        fftw_plan plan_backward_grady_delta = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(grady_delta_FT.data()), reinterpret_cast<fftw_complex*>(grady_delta.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_grady_delta);

	fftw_plan plan_backward_gradz_delta = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints, reinterpret_cast<fftw_complex*>(gradz_delta_FT.data()), reinterpret_cast<fftw_complex*>(gradz_delta.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_backward_gradz_delta);
	
        fftw_destroy_plan(plan_backward_gradx_delta);
        fftw_destroy_plan(plan_backward_grady_delta);
	fftw_destroy_plan(plan_backward_gradz_delta);
        fftw_cleanup();

	
      //compute integrand
        for (int i = 0; i < meshpoints; i++) {
            for (int j = 0; j < meshpoints; j++) {
                for (int k = 0; k < meshpoints; k++) {
                    // Compute the 1D index for (i, j, k)
                    int idx = i * meshpoints * meshpoints + j * meshpoints + k;
		    gradx_delta[idx] *= norm_factor;
		    grady_delta[idx] *= norm_factor;
		    gradz_delta[idx] *= norm_factor;

                    integrand[idx] = theta[idx] * delta[idx];

		    integrand[idx] += vx[idx] * gradx_delta[idx] + vy[idx] * grady_delta[idx] + vz[idx] * gradz_delta[idx];

                }
            }
        }
               
        fftw_plan forwardPlan = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints,
            reinterpret_cast<fftw_complex*>(integrand.data()),
            reinterpret_cast<fftw_complex*>(integrand.data()),
            FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(forwardPlan);
        fftw_destroy_plan(forwardPlan);

       
        for (int i = 0; i < meshpoints; i++) {
            for (int j = 0; j < meshpoints; j++) {
                for (int k = 0; k < meshpoints; k++) {
                    int index = i * meshpoints * meshpoints + j * meshpoints + k;
		    if (index > 0){
                    FTinverseLaplacian[index] = ( - integrand[index]* pow(boxsize, 2)) /
                        (pow(6.28318, 2) * (pow(std::min(i, meshpoints - i), 2) +
                            pow(std::min(j, meshpoints - j), 2) +
                            pow(std::min(k, meshpoints - k), 2) ));
		    }
                }
            }
        }
	//FTinverseLaplacian[0] = 0;
;
        // Create FFTW plan for inverse transform
        fftw_plan inversePlan = fftw_plan_dft_3d(meshpoints, meshpoints, meshpoints,
            reinterpret_cast<fftw_complex*>(FTinverseLaplacian.data()),
            reinterpret_cast<fftw_complex*>(integral.data()),
            FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(inversePlan);
        fftw_destroy_plan(inversePlan);

                
        for (int i = 0; i < meshpoints; i++) {
            for (int j = 0; j < meshpoints; j++) {
                for (int k = 0; k < meshpoints; ++k) {
                    int index = i * meshpoints * meshpoints + j * meshpoints + k;
                    integral[index] = integral[index].real() / (meshpoints * meshpoints * meshpoints);
                }
            }
        }
        for (int i = 0; i < meshpoints; i++) {
            for (int j = 0; j < meshpoints; j++) {
                for (int k = 0; k < meshpoints; k++) {
                    int idx = i * meshpoints * meshpoints + j * meshpoints + k;

                    piprimedatanext[idx] = -(Hubble / 2.0) * (gradsquared[idx]) * (timestep);

                    // Gradient squared terms
                    piprimedatanext[idx] += - gradpipiprime[idx] * timestep;

                    // Integral term
                    piprimedatanext[idx] += -(integral[idx] * timestep);

                }
            }
        }


	//Gaussian smoothing

	if (sigma > 0) {

	  fftw_plan plan_forward_piprimedatanext = fftw_plan_dft_3d(
								    meshpoints,
								    meshpoints,
								    meshpoints,
								    reinterpret_cast<fftw_complex*>(piprimedatanext.data()),
								    reinterpret_cast<fftw_complex*>(piprimedatanext_FT.data()),
								    FFTW_FORWARD,
								    FFTW_ESTIMATE
								    );

	  fftw_execute(plan_forward_piprimedatanext);

	  double k_cutoff = 6.28318/(sigma);

	  for (int i = 0; i < meshpoints; ++i) {
            for (int j = 0; j < meshpoints; ++j) {
	      for (int k = 0; k < meshpoints; ++k) {
		double k_mag = (6.28318 / boxsize) * std::sqrt(std::pow(std::min(i, meshpoints - i), 2) +
							       std::pow(std::min(j, meshpoints - j), 2) +
							       std::pow(std::min(k, meshpoints - k), 2));
		int idx = i * meshpoints * meshpoints + j * meshpoints + k;
		piprimedatanext_FT[idx] *= exp( - pow(k_mag, 2)/(2.0 * pow(k_cutoff, 2)) );
	      }
            }
	  }

	  fftw_plan plan_inverse_piprimedatanext = fftw_plan_dft_3d(
								    meshpoints,
								    meshpoints,
								    meshpoints,
								    reinterpret_cast<fftw_complex*>(piprimedatanext_FT.data()),
								    reinterpret_cast<fftw_complex*>(piprimedatanext.data()),
								    FFTW_BACKWARD,
								    FFTW_ESTIMATE
								    );
	  fftw_execute(plan_inverse_piprimedatanext);

	  for (int i = 0; i < meshpoints; ++i) {
            for (int j = 0; j < meshpoints; ++j) {
	      for (int k = 0; k < meshpoints; ++k) {
		int idx = i * meshpoints * meshpoints + j * meshpoints + k;
		piprimedatanext[idx] = piprimedatanext[idx].real() * norm_factor;
	      }
            }
	  }

        
	  // Clean up the FFTW plan
	  fftw_destroy_plan(plan_inverse_piprimedatanext);;
	  fftw_destroy_plan(plan_forward_piprimedatanext);
	  fftw_cleanup();

	} 

        // The 'integral' vector now contains the result of the inverse FFT

        for (int i = 0; i < meshpoints; i++) {
            for (int j = 0; j < meshpoints; j++) {
                for (int k = 0; k < meshpoints; k++) {
                    int idx = i * meshpoints * meshpoints + j * meshpoints + k;

                    pidatanext[idx] = (piprimedata[idx] * timestep);
                    piprimedatanext[idx] += ((-(2.0 * Hubble) * piprimedata[idx] - ( 2.0 * Hubble_prime ) * pidata[idx]) * timestep);

		    //EFT terms
		    piprimedatanext[idx] += (tildec2_s * theta[idx] * timestep);
		    piprimedatanext[idx] += (tildec2_v_overH * thetaprime[idx] * timestep);

		    pidata[idx] += pidatanext[idx];
		    piprimedata[idx] += piprimedatanext[idx];

                }
            }
        }

        time += timestep;
	cout << "current scale factor is: " << a_scale << "\n";
	lin_pi_growth_times_a2 += a_step/((H0/299792.458)*sqrt(((Omega_m/(pow(a_scale,3)))+(1-Omega_m))));
        a_scale += a_step;
        numsteps += 1;

        
    }
    cout << pidata[2*meshpoints*meshpoints + 3*meshpoints + 5];

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    
    cout <<endl<<endl<<"Running time for " << numsteps << " steps is :::: " << duration.count()<<" seconds.\n" << endl;

    cout << "Final time is " << time << " Mpc and final scale factor is " << a_scale << "\n";
   
     ofstream outputFile(pion_output_path);

    //Check if the file was opened successfully
    if (!outputFile) {
        cerr << "Error opening the file." << endl;
        return 1;
    }
    auto pidata2 = unflatten(meshpoints, pidata);
    outputFile << "[";
    for (int i = 0; i < meshpoints-1; i++) {
        outputFile << "[";
        for (int j = 0; j < meshpoints-1; j++) {
	  outputFile << "[";
            for (int k = 0; k < meshpoints - 1; k++) {
	      outputFile << pidata2[i][j][k].real() << ",";
            }
            outputFile << pidata2[i][j][meshpoints-1].real() << "],";
        }
	outputFile << "[";
        for (int k = 0; k < meshpoints - 1; k++) {
	   outputFile << pidata2[i][meshpoints-1][k].real() << ",";
            }
        outputFile << pidata2[i][meshpoints-1][meshpoints-1].real() << "]],";
    }
    outputFile << "[";
    for (int j = 0; j < meshpoints - 1; j++) {
       outputFile << "[";
            for (int k = 0; k < meshpoints - 1; k++) {
	      outputFile << pidata2[meshpoints - 1][j][k].real() << ",";
            }
            outputFile << pidata2[meshpoints - 1][j][meshpoints - 1].real() << "],";
        }
	outputFile << "[";
        for (int k = 0; k < meshpoints - 1; k++) {
	   outputFile << pidata2[meshpoints - 1][meshpoints-1][k].real() << ",";
            }
    outputFile << pidata2[meshpoints-1][meshpoints-1][meshpoints-1].real() << "]]]";

    // Close the file
    outputFile.close();

    ofstream outputFile2(pion_prime_output_path);

    //Check if the file was opened successfully
    if (!outputFile2) {
        cerr << "Error opening the file." << endl;
        return 1;
    }
    auto piprimedata2 = unflatten(meshpoints, piprimedata);
    outputFile2 << "[";
    for (int i = 0; i < meshpoints-1; i++) {
        outputFile2 << "[";
        for (int j = 0; j < meshpoints-1; j++) {
	  outputFile2 << "[";
            for (int k = 0; k < meshpoints - 1; k++) {
	      outputFile2 << piprimedata2[i][j][k].real() << ",";
            }
            outputFile2 << piprimedata2[i][j][meshpoints-1].real() << "],";
        }
	outputFile2 << "[";
        for (int k = 0; k < meshpoints - 1; k++) {
	   outputFile2 << piprimedata2[i][meshpoints-1][k].real() << ",";
            }
        outputFile2 << piprimedata2[i][meshpoints-1][meshpoints-1].real() << "]],";
    }
    outputFile2 << "[";
    for (int j = 0; j < meshpoints - 1; j++) {
       outputFile2 << "[";
            for (int k = 0; k < meshpoints - 1; k++) {
	      outputFile2 << piprimedata2[meshpoints - 1][j][k].real() << ",";
            }
            outputFile2 << piprimedata2[meshpoints - 1][j][meshpoints - 1].real() << "],";
        }
	outputFile2 << "[";
        for (int k = 0; k < meshpoints - 1; k++) {
	   outputFile2 << piprimedata2[meshpoints - 1][meshpoints-1][k].real() << ",";
            }
    outputFile2 << piprimedata2[meshpoints-1][meshpoints-1][meshpoints-1].real() << "]]]";

    // Close the file
    outputFile2.close();

    ofstream outputFile3(overdensity_output_path);

    //Check if the file was opened successfully
    if (!outputFile3) {
        cerr << "Error opening the file." << endl;
        return 1;
    }
    auto delta2 = unflatten(meshpoints, overdensity);
    outputFile3 << "[";
    for (int i = 0; i < meshpoints-1; i++) {
        outputFile3 << "[";
        for (int j = 0; j < meshpoints-1; j++) {
	  outputFile3 << "[";
            for (int k = 0; k < meshpoints - 1; k++) {
	      outputFile3 << delta2[i][j][k].real() << ",";
            }
            outputFile3 << delta2[i][j][meshpoints-1].real() << "],";
        }
	outputFile3 << "[";
        for (int k = 0; k < meshpoints - 1; k++) {
	   outputFile3 << delta2[i][meshpoints-1][k].real() << ",";
            }
        outputFile3 << delta2[i][meshpoints-1][meshpoints-1].real() << "]],";
    }
    outputFile3 << "[";
    for (int j = 0; j < meshpoints - 1; j++) {
       outputFile3 << "[";
            for (int k = 0; k < meshpoints - 1; k++) {
	      outputFile3 << delta2[meshpoints - 1][j][k].real() << ",";
            }
            outputFile3 << delta2[meshpoints - 1][j][meshpoints - 1].real() << "],";
        }
	outputFile3 << "[";
        for (int k = 0; k < meshpoints - 1; k++) {
	   outputFile3 << delta2[meshpoints - 1][meshpoints-1][k].real() << ",";
            }
    outputFile3 << delta2[meshpoints-1][meshpoints-1][meshpoints-1].real() << "]]]";

    // Close the file
    outputFile3.close();

    cout << "Output files written: " << pion_output_path << " " << pion_prime_output_path << " " << overdensity_output_path << "\n";

    

    return 0;
}

