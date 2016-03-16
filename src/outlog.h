#include<iostream>
#include<fstream>
using namespace std;

class outlog:public ostream {
        private:
            std::ofstream logfile;
        public:
            outlog(const char* file)
            {
                logfile.open(file);

                if (!(logfile.is_open())) {
                    std::cerr << "[]: Couldn't open file \"" << file << "\" for logging.\n";
                    this->~outlog(); /* Destroy the object */
                }
            }

            ~outlog()
            {
                if (logfile.is_open())
                    logfile.close();
            }

            template <class T>
            outlog& operator<<(const T& out)
            {
                std::cout << out;
                logfile << out;

                return *this;
            }
	streamsize  precision(streamsize prec)
	{
		std::cout.precision(prec);
		return logfile.precision(prec);
		
	}

        // this is the type of std::cout
    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;

    // this is the function signature of std::endl
    typedef CoutType& (*StandardEndLine)(CoutType&);

    // define an operator<< to take in std::endl
    outlog& operator<<(StandardEndLine manip)
    {
        // call the function, but we cannot return it's value
        manip(std::cout);
        manip(logfile);

        return *this;
    }
    
    };


extern outlog log_cout;
#define cout log_cout





