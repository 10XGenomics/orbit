#ifndef PARAMETERSINFO_DEF
#define PARAMETERSINFO_DEF

class ParameterInfoBase {
public:
    const char* nameString; //string that identifies parameter
    int inputLevel; //where the parameter was defined
    int inputLevelAllowed; //at which inpurt level parameter definition is allowed
    virtual void inputValues(istringstream &streamIn) =0;
    friend std::ostream& operator<< (std::ostream& o, ParameterInfoBase const& b);
    ParameterInfoBase(const char* nameString, int inputLevel, int inputLevelAllowed);
    virtual ~ParameterInfoBase() = default;
    ParameterInfoBase(const ParameterInfoBase&) = default;
    ParameterInfoBase &operator=(const ParameterInfoBase&) = delete;
    ParameterInfoBase(ParameterInfoBase&&) = default;	
    ParameterInfoBase& operator=(ParameterInfoBase&&) = default;
protected:
    virtual void printValues(std::ostream& o) const = 0;
};



inline std::ostream& operator<< (std::ostream& o, ParameterInfoBase const& b) {
    b.printValues(o);
    return o;
};


template <class parameterType>
inline parameterType inputOneValue (istringstream &streamIn) {
    parameterType oneV;
    streamIn >> oneV;
    return oneV;
};
template <>
inline string inputOneValue <string> (istringstream &streamIn) {
    string oneV="";
    streamIn >> ws;//skip whitespace
    if (streamIn.peek()!='"') {//simple parameter with no spaces or "
        streamIn >> oneV;
    } else {
        streamIn.get();//skip "
        getline(streamIn,oneV,'"');
    };
    return oneV;
};
// Make sure -1 gets wrapped to UINT_MAX.  It is nessessary to do it this way
// in case of older versions of libc.
template <>
inline uint inputOneValue <uint> (istringstream &streamIn) {
    int64_t oneV;
    streamIn >> oneV;
    return static_cast<uint>(oneV);
};


template <class parameterType>
inline void printOneValue (parameterType *value, std::ostream& outStr) {
    outStr << *value;
};
template <>
inline void printOneValue <string> (string *value, std::ostream& outStr) {
    if ((*value).find_first_of(" \t")!=std::string::npos) {//there is white space in the argument, put "" around
        outStr << '\"' << *value <<'\"';
    } else {
        outStr << *value;
    };
};

template <class parameterType>
class ParameterInfoScalar : public ParameterInfoBase {
public:
    parameterType *value;
    vector <parameterType> allowedValues;

    ParameterInfoScalar(int inputLevelIn, int inputLevelAllowedIn, const char* nameStringIn, parameterType* valueIn)
        : ParameterInfoBase(nameStringIn, inputLevelIn, inputLevelAllowedIn),
          value(valueIn) {};

    void inputValues(istringstream &streamIn) override {
        *value=inputOneValue <parameterType> (streamIn);
    };

protected:
   void printValues(std::ostream& outStr) const override {
       printOneValue(value, outStr);
   };

};

template <class parameterType>
class ParameterInfoVector : public ParameterInfoBase {
public:
    vector <parameterType> *value;
    vector <parameterType> allowedValues;

    ParameterInfoVector(int inputLevelIn, int inputLevelAllowedIn, const char* nameStringIn, vector <parameterType> *valueIn)
        : ParameterInfoBase(nameStringIn, inputLevelIn, inputLevelAllowedIn),
          value(valueIn) {};

    void inputValues(istringstream &streamIn) override {
        (*value).clear();
        while (streamIn.good()) {
            (*value).push_back(inputOneValue <parameterType> (streamIn));
            streamIn >> ws; //remove white space, may arrive at the end of line
        };
    };

protected:
   void printValues(std::ostream& outStr) const override {
       for (int ii=0; ii < (int) (*value).size(); ii++) {
           printOneValue(&(*value).at(ii),outStr);
           outStr<<"   ";
       };
   };
};
#endif
