#include <iostream>
#include <amp.h>
#include <amp_math.h>
#include <vector>
#include <functional>
using namespace Concurrency;
class A
{ 
public:
    A(const array_view<int>& arr) : a(arr) {}
    void operator() (index<1> idx) restrict(amp)
    {
        a[idx] += (1);
    }

protected:
    array_view<int> a;
};

int main()
{
    std::vector<int> arr = { 1,2,3,4,5 };
    array_view<int> arr_view(arr.size(), arr);
    parallel_for_each(arr_view.extent,
        A(arr_view));
    arr_view.synchronize();
    for (int i = 0; i < arr.size(); ++i)
    {
        std::cout << arr[i] << std::endl;
    }
    std::getchar();
    return 0;
}