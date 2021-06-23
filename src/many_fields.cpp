#include <iostream>
#include <H5Cpp.h>
#include "field.hpp"


#ifndef _H5_NO_NAMESPACE_
using namespace H5;
#ifndef _H5_NO_STD_
    using std::cout;
    using std::endl;
#endif /* _H5_NO_STD_ */        
#endif /* _H5_NO_NAMESPACE_ */



System::System() : Field()
{
    unsigned seed;
    seed = time(0);
    generator.seed(seed);

    float x0s[5] = {-40, -20, 0, 20, 40};
    float y0s[3] = {-60, -40, -20};
    int index = 0;
    num_of_field = 15;
    sys = new Field [num_of_field];

    double v0[2], angle;
    for(int i=0; i<5; i++)
        for(int j=0; j<3; j++)
        {   
            std::uniform_real_distribution<double> uniform(-PI, PI);
            angle = uniform(generator);
            float vv = 0.035;
            v0[0] = 0.035 * cos(angle);
            v0[1] = 0.035 * sin(angle);
            sys[index] = Field(x0s[i], y0s[j], 12, v0);  

            unsigned seed;
            seed = time(0) * (index+1);
            sys[index].generator.seed(seed);
            index += 1;   
        }

    confinement = dmatrix(0, fullN_0, 0, fullN_1);
    for(int i=70; i<=fullN_0; i++)
        for(int j=0; j<=40; j++)
            confinement[i][j] = 1.;

    for(int i=70; i<=fullN_0; i++)
        for(int j=60; j<=fullN_1; j++)
            confinement[i][j] = 1.;

    for(int i=0; i<=10; i++)
        for(int j=0; j<=fullN_1; j++)
            confinement[i][j] = 1.;

    conf = dmatrix(0, fullN_0, 0, fullN_1);
    for(int i=70; i<=fullN_0; i++)
        for(int j=0; j<=fullN_1; j++)
            conf[i][j] = 1.;

    for(int i=0; i<=10; i++)
        for(int j=0; j<=fullN_1; j++)
            conf[i][j] = 1.;
}

void System::save_data(const char *str)
{
    char fileName[50];
    sprintf(fileName, "%s.csv", str);

    std::ifstream fin(fileName);
    if (!(fin))
    {
        std::ofstream fout(fileName);
        
        double ** phi;
        phi = dmatrix(0, fullN_0, 0, fullN_1);
        for (int i = 0; i <= fullN_0; i++)
            for (int j = 0; j <= fullN_1; j++)
            {
                phi[i][j] = 0.;
                for (int k=0; k<num_of_field; k++)
                    phi[i][j] += sys[k].fullphi[i][j];

                phi[i][j] += confinement[i][j];
            }
            
                
        for (int i = 0; i <= fullN_0; i++)
        {
            for (int j = 0; j <= fullN_1; j++)
            {
                fout << phi[i][j];

                if (j != fullN_1)
                {
                    fout << ",";
                }
                else
                {
                    fout << "\n";
                }
            }
        }
        free_dmatrix(phi, 0, fullN_0, 0, fullN_1);
        fout.close();
    }
}

void System::update(float dt, bool therm, float k_theta, float omega)
{
    double **h, **h1, **h1_laplace, **temp;
    h = dmatrix(0, fullN_0, 0, fullN_1);
    h1 = dmatrix(0, fullN_0, 0, fullN_1);
    h1_laplace = dmatrix(0, fullN_0, 0, fullN_1);
    temp = dmatrix(0, fullN_0, 0, fullN_1);

    for (int i = 0; i <= fullN_0; i++)
        for (int j = 0; j <= fullN_1; j++)
        {
            h[i][j] = 0.;
            h1[i][j] = 0.;
            h1_laplace[i][j] = 0.;
            temp[i][j] = 0.;
        }

    for(int k=0; k<num_of_field; k++)
        for (int i = 0; i <= fullN_0; i++)
            for (int j = 0; j <= fullN_1; j++)
            {
                double val = sys[k].fullphi[i][j];
                h[i][j] += val * val;
                h1[i][j] += val;
            }

    second_diff(fullN_0, fullN_1, h1, 0, h1_laplace);
    second_diff(fullN_0, fullN_1, h1, 1, temp);
    for (int i = 0; i <= fullN_0; i++)
        for (int j = 0; j <= fullN_1; j++)
            h1_laplace[i][j] += temp[i][j];

    if(therm == false)
    {
        for(int k=0; k<num_of_field; k++)
            sys[k].update(dt, h, h1_laplace, confinement,
                    true, true, true, true, k_theta, omega);
    }
    else
    {
        for(int k=0; k<num_of_field; k++)
            sys[k].update(dt, h, h1_laplace, conf,
                    true, true, true, false, k_theta, omega);
    }
    
    free_dmatrix(h, 0, fullN_0, 0, fullN_1);
    free_dmatrix(h1, 0, fullN_0, 0, fullN_1);
    free_dmatrix(h1_laplace, 0, fullN_0, 0, fullN_1);
    free_dmatrix(temp, 0, fullN_0, 0, fullN_1);
    
}


int System::simulation(const char *str, float T, float dt, bool therm, float k_theta, float omega)
{   
    float t, ti(0);
    t = ti;
    int k=0;

    char fileName[50];
    sprintf(fileName, "%s.h5", str);

    const H5std_string FILE_NAME(fileName);
    const H5std_string ATTR_NAME1("subdomain");
    const H5std_string ATTR_NAME2("center");
    const H5std_string ATTR_NAME3("velocity"); 
    try
    {
        /* 
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately.
         */
        Exception::dontPrint();

        // Create a new file using the default property lists.
        // H5::H5F_ACC_TRUNC : create a new file or overwrite an existing file.
        H5File file(FILE_NAME, H5F_ACC_TRUNC); 

        const H5std_string DATASET_NAME("confinement");
        // Use H5::hsize_t (similar to int) for dimensions.
        hsize_t dims[2];               // dataset dimensions
        dims[0] = fullN_0+1;
        dims[1] = fullN_1+1;
        
        double data[fullN_0+1][fullN_1+1];
        if(therm == false)
        {
            for(int i=0; i<=fullN_0; i++)
                for(int j=0; j<=fullN_1; j++)
                    data[i][j] = confinement[i][j];
        }
        else
        {  
            for(int i=0; i<=fullN_0; i++)
                for(int j=0; j<=fullN_1; j++)
                    data[i][j] = conf[i][j];
        }
        // Create the dataspace for a dataset first.
        DataSpace dataspace(2, dims);
        
        // Create the dataset under group with specified dataspace.      
        DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);

        // Write data in buffer to dataset.
        dataset.write(data, PredType::NATIVE_DOUBLE);

        for(int kk=0; kk<num_of_field; kk++)
        {   
            char groupName[50];
            sprintf(groupName, "field_%d", kk);


            char datasetName[50];
            sprintf(datasetName, "t_%g", t);

            const H5std_string GROUP_NAME(groupName);
            const H5std_string DATASET_NAME(datasetName);

            // Create a group under root '/'.
            Group group(file.createGroup(GROUP_NAME));
            
            const int RANK = 2;
            // Use H5::hsize_t (similar to int) for dimensions.
            hsize_t dims[RANK];               // dataset dimensions
            dims[0] = N+1;
            dims[1] = N+1;
            
            double data[N+1][N+1];
            for(int i=0; i<=N; i++)
                for(int j=0; j<=N; j++)
                    data[i][j] = sys[kk].phi[i][j];

            // Create the dataspace for a dataset first.
            DataSpace dataspace(RANK, dims);
            
            // Create the dataset under group with specified dataspace.      
            DataSet dataset = group.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);

            // Write data in buffer to dataset.
            dataset.write(data, PredType::NATIVE_DOUBLE);

            /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
            hsize_t attr1_dims[1] = {2};             // attribute dimension, rank = 1
            
            // Create the dataspace for an attribute first.
            DataSpace attr1_dataspace(1, attr1_dims);         // rank = 1

            // Create the attribute of dataset with specified dataspace.
            Attribute attribute1 = dataset.createAttribute(ATTR_NAME1, PredType::STD_I32BE, attr1_dataspace);
            
            // Write data in buffer to attribute.
            attribute1.write(PredType::NATIVE_INT, sys[kk].subdomain);


            hsize_t attr2_dims[1] = {2};
            DataSpace attr2_dataspace(1, attr2_dims);
            Attribute attribute2 = dataset.createAttribute(ATTR_NAME2, PredType::NATIVE_DOUBLE, attr2_dataspace);
            attribute2.write(PredType::NATIVE_DOUBLE, sys[kk].center);

            hsize_t attr3_dims[1] = {2};
            DataSpace attr3_dataspace(1, attr3_dims);
            Attribute attribute3 = dataset.createAttribute(ATTR_NAME3, PredType::NATIVE_DOUBLE, attr3_dataspace);
            attribute3.write(PredType::NATIVE_DOUBLE, sys[kk].v0);
            
            // Save and exit the group.
            group.close();
        }
        // Save and exit the file.
        file.close();

    }  // end of try block

    
    // Catch failure caused by the H5File operations.
    catch(FileIException error)
    {
        error.printErrorStack();
        return -1;
    }

    // Catch failure caused by the DataSet operations.
    catch(DataSetIException error)
    {
        error.printErrorStack();
        return -1;
    }

    // Catch failure caused by the DataSpace operations.
    catch(DataSpaceIException error)
    {
        error.printErrorStack();
        return -1;
    }

    while(t <= T+ti)
    {
        k += 1;
        t += dt;
        update(dt, therm, k_theta, omega);
        
        if(k % 20 == 0)
        {

            try
            {
                /* 
                 * Turn off the auto-printing when failure occurs so that we can
                 * handle the errors appropriately.
                 */
                Exception::dontPrint();

                H5File file(FILE_NAME, H5F_ACC_RDWR); 

                for(int kk=0; kk<num_of_field; kk++)
                {   
                    char groupName[50];
                    sprintf(groupName, "field_%d", kk);


                    char datasetName[50];
                    sprintf(datasetName, "t_%g", t);

                    const H5std_string GROUP_NAME(groupName);
                    const H5std_string DATASET_NAME(datasetName);

                    // Create a group under root '/'.
                    Group group = file.openGroup(GROUP_NAME);
                    
                    const int RANK = 2;
                    // Use H5::hsize_t (similar to int) for dimensions.
                    hsize_t dims[RANK];               // dataset dimensions
                    dims[0] = N+1;
                    dims[1] = N+1;
                    
                    double data[N+1][N+1];
                    for(int i=0; i<=N; i++)
                        for(int j=0; j<=N; j++)
                            data[i][j] = sys[kk].phi[i][j];

                    // Create the dataspace for a dataset first.
                    DataSpace dataspace(RANK, dims);
                    
                    // Create the dataset under group with specified dataspace.      
                    DataSet dataset = group.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);

                    // Write data in buffer to dataset.
                    dataset.write(data, PredType::NATIVE_DOUBLE);

                    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
                    hsize_t attr1_dims[1] = {2};             // attribute dimension, rank = 1
                    
                    // Create the dataspace for an attribute first.
                    DataSpace attr1_dataspace(1, attr1_dims);         // rank = 1

                    // Create the attribute of dataset with specified dataspace.
                    Attribute attribute1 = dataset.createAttribute(ATTR_NAME1, PredType::STD_I32BE, attr1_dataspace);
                    
                    // Write data in buffer to attribute.
                    attribute1.write(PredType::NATIVE_INT, sys[kk].subdomain);

       
                    hsize_t attr2_dims[1] = {2};
                    DataSpace attr2_dataspace(1, attr2_dims);
                    Attribute attribute2 = dataset.createAttribute(ATTR_NAME2, PredType::NATIVE_DOUBLE, attr2_dataspace);
                    attribute2.write(PredType::NATIVE_DOUBLE, sys[kk].center);
       
                    hsize_t attr3_dims[1] = {2};
                    DataSpace attr3_dataspace(1, attr3_dims);
                    Attribute attribute3 = dataset.createAttribute(ATTR_NAME3, PredType::NATIVE_DOUBLE, attr3_dataspace);
                    attribute3.write(PredType::NATIVE_DOUBLE, sys[kk].v0);
                    
                    // Save and exit the group.
                    group.close();
                }
                // Save and exit the file.
                file.close();

            }  // end of try block

            
            // Catch failure caused by the H5File operations.
            catch(FileIException error)
            {
                error.printErrorStack();
                return -1;
            }

            // Catch failure caused by the DataSet operations.
            catch(DataSetIException error)
            {
                error.printErrorStack();
                return -1;
            }

            // Catch failure caused by the DataSpace operations.
            catch(DataSpaceIException error)
            {
                error.printErrorStack();
                return -1;
            }
		/*
        if(therm == true)
            procBar(int(k * dt * 100/T), "Thermalize");
        else
            procBar(int(k * dt * 100/T));
		*/
        }
    }
    //std::cout<<std::endl;
    return 0;
}

