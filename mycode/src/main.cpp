#include "admm.h"
#include "other/properties.h"

int main(int argc, char **argv) {
    Properties properties(argc, argv);
    int dimension = properties.GetInt("dimension");
    int min_barrier = properties.GetInt("min_barrier");
    int max_delay = properties.GetInt("max_delay");
    int max_iterations = properties.GetInt("max_iterations");
    int interval = properties.GetInt("interval");
	int filter_type= properties.GetInt("filter_type");
    double rho = properties.GetDouble("rho");
    double l2reg = properties.GetDouble("l2reg");
	double l1reg = properties.GetDouble("l1reg");
    double ABSTOL = properties.GetDouble("ABSTOL");
    double RELTOL = properties.GetDouble("RELTOL");
	double theta= properties.GetDouble("theta");
	
    std::string train_data_path = properties.GetString("train_data_path");
    std::string test_data_path = properties.GetString("test_data_path");
	
	if(argc>=2)
	{
		filter_type=atoi(argv[1]);
		min_barrier=atoi(argv[2]);
		//std::cout<<filter_type<<"##";
	}
	
	 
    ADMM admm(dimension, min_barrier, max_delay, max_iterations, interval, rho, l2reg, l1reg, ABSTOL, RELTOL, train_data_path,
              test_data_path,filter_type,theta);
    if (admm.GetID() == 0) {
		std::cout<<filter_type<<"**"<<min_barrier<<std::endl;
        properties.Print();
    }
    admm.Run();
    return 0;
}
