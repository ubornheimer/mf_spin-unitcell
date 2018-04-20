#include "../include/Parameters.h"

Parameters::Parameters()
{}

Parameters::Parameters(int argc, char* argv[])
{
    if (argc != 12)
    {
        printf("command-line error!! argc = %i\n", argc);
        exit(1);
    }

    unitcellnr = atoi(argv[1]);
    linknr = atoi(argv[2]);
    compnr = atoi(argv[3]);
    repnr = atoi(argv[4]);
    gentype = atoi(argv[5]);
    a0 = atof(argv[6]);
    a2 = atof(argv[7]);
    tol = atoi(argv[8]);
    maxit = atoi(argv[9]);
    minit = atoi(argv[10]);
    sitenr = atoi(argv[11]);

    phi_in = read_phiin();
    coeffgen = read_coeffgen();
}

Parameters::~Parameters()
{}

cd3D Parameters::read_coeffgen()
{
    cd3D result;

    for (int i=0; i<linknr; i++)
    {
        string filename = "coeffgen"+std::to_string(i)+".txt";
        result.push_back(ffile_coeffgen(filename));
    }

    return result;
}

cd2D Parameters::ffile_coeffgen(string filename)
{
    cd2D result;
    cd1D temp;

    ifstream input;
    input.open(filename);
    bool isdefault = false;
    int counter = 0;

    if(input.good())
    {
        string inputline;

        while(getline(input,inputline))
        {
            string buffer;
            stringstream ss(inputline);
            vector<string> values;

            while(ss >> buffer)
                values.push_back(buffer);

            if(values.size() == 2)
            {
                //printf("value size good!");
                if (counter%repnr==0 and counter>0)
                {
                    result.push_back(temp);
                    temp.clear();
                    counter = 0;
                    //printf("\n");
                }
                cd test = (cd)atof(values[0].c_str())/a0 + (cd)atof(values[1].c_str())/a2;
                //printf("%f\n", real(test));
                temp.push_back(test);
                counter++;
            }
            else
            {
                isdefault = true;
                //printf("value size not good");
            }
        }
        result.push_back(temp);
    }
    else
    {
        isdefault = true;
        //printf("input not good");
    }
    input.close();

    if(isdefault)
    {
        printf("Initial read of the coefficients of the generators from file failed. Default parameters used.\n");
        for (int i=0; i<repnr; i++)
        {
            if (i>0)
            {
                result.push_back(temp);
                temp.clear();
            }
            for(int j=0; j<repnr; j++)
                temp.push_back((cd)0.05);
        }
        result.push_back(temp);
    }
    else
        printf("Initial read of the generator coefficients from file succeeded.\n");

    return result;
}

cd4D Parameters::read_phiin() //csize = cluster size, repnr = nr of elements in representation
{
    cd4D result;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0,1);

    for (int i=0; i<sitenr; i++)
    {
        cd3D ttttemp;
        for (int j=0; j<sitenr; j++)
        {
            cd2D tttemp;
            for (int k=0; k<unitcellnr; k++)
            {
                cd1D ttemp;
                for(int j=0; j<repnr; j++)
                {
                    cd temp = dis(gen);
                    ttemp.push_back(temp);
//                    printf("%f \t",real(temp));
                }
                tttemp.push_back(ttemp);
//                printf("\n");
            }
            ttttemp.push_back(tttemp);
//            printf("-----------------\n");
        }
        result.push_back(ttttemp);
//        printf("-----------------\n");
//        printf("-----------------\n");
    }

    return result;
}
