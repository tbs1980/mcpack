/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_TEXTDATAIO_HPP
#define MCPACK_TEXTDATAIO_HPP

namespace mcpack{ namespace utils {


    //an exception class for file opening error
    class TextDataException: public std::exception
    {
    public:
        TextDataException()throw ()
        {}

        TextDataException(std::string message) throw () 
        :m_message(message) 
        {}

        virtual ~TextDataException()throw () 
        {} 

        virtual const char* what() const throw()
        {
            return m_message.c_str();
        }

    private:
        std::string m_message;
    };


    //conversion from string to number
    //see http://www.cplusplus.com/articles/D9j2Nwbp/

    template <typename T>
    T StringToNumber (  std::string const & Text )
    {
        std::istringstream ss(Text);
        T result;
        return ss >> result ? result : 0;
    }


    
    /*
    This function is the base for reading a matrix<T> from a file. 
    The entries are supposed to be tabe or space separated. 
    All the rows should have the same number of entries (columns).
    */
    
    template <typename T>
    void ReadMatrixFromTextFile(std::vector< std::vector<T> > & cont,std::string fileName,std::string separation)
    {
        //try to open file
        std::ifstream inFile;
        inFile.open(fileName.c_str(),std::ios::in);

        unsigned lineNo=0;
        unsigned ncols=0;
        if(inFile.is_open())
        {
            //clear the contents first
            //cont.clear();

            while(!inFile.eof())
            {
                std::string line;
                std::getline(inFile,line);
                std::istringstream ss(line);
                std::string entry;
                std::vector<T> row;

                while(std::getline(ss,entry,*(separation.c_str()) ) )
                {
                    T num=StringToNumber<T>(entry);
                    row.push_back(num);
                }
                    //check if all the rows have the same number of columns
                if(lineNo==0)
                {
                    ncols=row.size();
                    cont.push_back(row);
                }
                else
                {
                    //if all rows does not have the same number of columns 
                    //stop reading. The last line will not have any elements.

                    if(row.size()==ncols) 
                    {
                        cont.push_back(row);
                    }
                    else
                    {
                        if(row.size()==0) 
                        {
                            break;// we are most likely at the end of the file
                        }
                        else
                        {                           
                            std::string message=std::string("\nIn file ")+fileName+std::string(", all rows do not have the same number of entries  ");
                            throw TextDataException(message);
                        }
                    }
                }               
                ++lineNo;
            }           
                inFile.close();     
        }
        else
        {
            std::string message=std::string("Error in opening the file ")+fileName;
            throw TextDataException(message);
        }
    }

    //same as the above function, but for a case where no delimiter is present
    template <typename T>
    void ReadMatrixFromTextFile(std::vector< std::vector<T> > & cont,std::string fileName)
    {
        //try to open file
        std::ifstream inFile;
        inFile.open(fileName.c_str(),std::ios::in);

        unsigned lineNo=0;
        unsigned ncols=0;
        if(inFile.is_open())
        {
            //clear the contents first
            //cont.clear();

            while(!inFile.eof())
            {
                std::string line;
                std::getline(inFile,line);
                std::istringstream ss(line);
                //std::string entry;
                std::vector<T> row;

                while(!ss.eof())
                {
                    //ss>>entry;
                    //T num=StringToNumber<T>(entry);
                    T num;
                    ss>>num;
                    row.push_back(num);
                }

                //check if all the rows have the same number of columns
                if(lineNo==0)
                {
                    ncols=row.size();
                    cont.push_back(row);
                }
                else
                {
                    //if all rows does not have the same number of columns 
                    //stop reading. The last line will not have any elements.

                    if(row.size()==ncols) 
                    {
                        cont.push_back(row);
                    }
                    else
                    {
                        if(row.size()==0) 
                        {
                            break;// we are most likely at the end of the file
                        }
                        else
                        {   
                            //this is slightly different from the earlier version
                            //this is due to the fact that we are using the oprator >>
                            //instead of getline()
                            if(inFile.eof())
                            {
                                break;
                            }
                            else
                            {               
                                std::string message=std::string("In file ")+fileName+std::string(", all rows do not have the same number of entries  ");
                                throw TextDataException(message);
                            }
                        }
                    }
                }               
                ++lineNo;
            }           
                inFile.close();     
        }
        else
        {
            std::string message=std::string("Error in opening the file ")+fileName;
            throw TextDataException(message);
        }
    }





    //A function to write matrix into file
    //only works with Eigen as we use rows() and cols() functions
    template<class realMatrixType>
    void WriteMatrix2TextFile(realMatrixType const & mat,std::string fileName,
        unsigned precision,std::string separation)
    {
        std::ofstream file;
        file.open(fileName.c_str(),std::ios::trunc);

        if(file.is_open())
        {
   
            file<<std::scientific;
            for(unsigned i=0;i<mat.rows();++i)
            {
                for(unsigned j=0;j<mat.cols()-1;++j)
                {
                    file<<std::setprecision(precision)<<mat(i,j)<<separation;
                }
                file<<std::setprecision(precision)<<mat(i,(mat.cols()-1) )<<std::endl;
            }
            file.close();
        }
        else
        {
            std::string message=std::string("Error in opening the file ")+fileName;
            throw TextDataException(message);
        }
    }


    
    //A function to read matrix from a text file
    //only works with Eigen as we use rows() and cols() functions
    template<class realMatrixType>
    void ReadMatrixFromTextFile(realMatrixType & mat,std::string fileName,std::string separation)
    {
        typedef typename realMatrixType::Scalar scalarType;
        //first read matrix into container
        std::vector< std::vector<scalarType> > container;
        typedef typename std::vector< std::vector<scalarType> >::size_type sizeType;
        ReadMatrixFromTextFile<scalarType>(container,fileName,separation);
            
        //now check if the sizes are the same
        if(container.size()!= (sizeType)mat.rows() || container[0].size()!= (sizeType) mat.cols())
        {
            std::cout<<container.size()<<"\t"<< (sizeType)mat.rows()<<"\t"<<container[0].size() <<"\t"<<(sizeType) mat.cols()<<std::endl;
            std::string message=std::string("Dimensions of the matrix do not agree with the dimension of the matrix in ")+fileName;
            throw TextDataException(message);           
        }
        else
        {
            //assign values
            for(int i=0;i<mat.rows();++i)
            {
                for(int j=0;j<mat.cols();++j)
                {
                    mat(i,j)=container[i][j];
                }
            }           
        }
            
    }

    template<class RealVectorType>
    RealVectorType String2Vector(std::string str,std::string separation)
    {
        typedef typename RealVectorType::Index IndexType;
        typedef typename RealVectorType::Scalar RealType;
        std::istringstream ss(str);
        std::string entry;
        std::vector<RealType> row;

        while(std::getline(ss,entry,*(separation.c_str()) ) )
        {
            RealType num=mcpack::utils::StringToNumber<RealType>(entry);
            row.push_back(num);
        }

        RealVectorType vect(row.size());
        for(IndexType i=0;i<(IndexType)row.size();++i)
        {
            vect(i)=row[i];
        }
        
        return vect;
    }
    
} //namespace utils

} //namespace mcpack

#endif //MCPACK_TEXTDATAIO_HPP
