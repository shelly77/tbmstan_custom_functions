
//C++ template to define the statistical model.


#include <TMB.hpp>

template <class Type> //convert x into range (0,1)
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);} 


 template<class Type> // define generalised poisson distribution
 inline Type dgpois(const Type &x, const Type &lambda, const Type &omega, int give_log=0)
 {
  Type lambdastar = (Type(1) - omega)*lambda + omega*x;
  Type logres = log(1 - omega) + log(lambda) + (x - 1)*log(lambdastar) -
     lgamma(x + 1) - lambdastar;
   if (give_log) return logres; else return exp(logres);

 }


template<class Type> // define zero-inflated binomial distribution
inline Type dzibinom(const Type &x, const Type &size, const Type &p, const Type & zip,
                     int give_log=0)
{
  Type logres;
  if (x==Type(0)) logres=log(zip + (Type(1)-zip)*dbinom(x, size, p, false)); 
  else logres=log(Type(1)-zip) + dbinom(x, size, p, true);
  if (give_log) return logres; else return exp(logres);
}

template<class Type> // define beta-binomial distribution
inline Type dbetabinom(const Type &x, const Type &size, const Type &alpha, const Type &beta, 
                       int give_log=0) // size is the total number of trials and x is the number of success
{
  Type logres =
    lfactorial(size)            
  - lfactorial(x)               
  - lfactorial(size - x)        
  + lgamma(alpha + x)  
  + lgamma(alpha + beta)
  - lgamma(alpha + beta + size)
  - lgamma(alpha)  
  - lgamma(beta)
  + lgamma(beta + size - x+ 1e-100);  
  if (give_log) return logres; else return exp(logres);
}

template<class Type> // define zero-inflated beta-binomial distribution
inline Type dzibetabinom(const Type &x, const Type &size, const Type &alpha, const Type &beta, const Type & zip,
                         int give_log=0) // size is the total number of trials and x is the number of success
{
  Type logres;
  if (x==Type(0)) logres=log(zip + (Type(1)-zip)*dbetabinom(x, size, alpha, beta, false)); 
  else logres=log(Type(1)-zip) + dbetabinom(x, size, alpha, beta, true);
  if (give_log) return logres; else return exp(logres);
}

