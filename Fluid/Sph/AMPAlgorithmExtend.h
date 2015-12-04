#pragma once
#ifndef SPH_AMP_ALGORITHM_EXTEND_H
#define SPH_AMP_ALGORITHM_EXTEND_H
#include <amp.h>
namespace Fluid
{
    template<typename DataType>
    extern void scan(Concurrency::array<DataType>& arr)
    {
        Concurrency::array<DataType> output(arr.extent, arr.accelerator_view);
        int size = arr.extent[0];
        int step = 1;
        while (step < size)
        {
            parallel_for_each(arr.extent,
                [&arr, &output, step](Concurrency::index<1> idx) restrict(amp)
            {
                int i = idx[0];
                if (i - step >= 0)
                {
                    output[i] = arr[i] + arr[i - step];
                }
                else
                {
                    output[i] = arr[i];
                }
            });
            parallel_for_each(arr.extent,
                [&arr, &output](Concurrency::index<1> idx) restrict(amp)
            {
                arr[idx] = output[idx];
            });
            step *= 2;
        }
    }

    template <typename DataType, typename FunctionType>
    extern DataType reduce(Concurrency::array<DataType>& arr, FunctionType function)
    {
        Concurrency::array<DataType> result(arr.extent, arr.accelerator_view);
        Concurrency::array<DataType> tmp(arr.extent, arr.accelerator_view);
        arr.copy_to(tmp);
        int size = arr.extent.size();
        while (size > 1)
        {
            int half_size = (size - 1) / 2 + 1;
            Concurrency::parallel_for_each(Concurrency::extent<1>(half_size),
                [size, &result](index<1> idx) restrict(amp)
            {
                result[idx] = tmp[idx * 2];
                if ((idx * 2 + 1)[0] < size)
                    result[idx] += tmp[idx * 2 + 1];
            });
            size = half_size;
            result.copy_to(tmp);
        }
    }
}
#endif