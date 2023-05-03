#include "infrastructure/RS_polynomial.h"
#include <thread>
#include <iostream>
prime_field::field_element* __dst[3];
prime_field::field_element* twiddle_factor;

void init_scratch_pad(int order)
{
    __dst[0] = new prime_field::field_element[order];//(prime_field::field_element*)malloc(order * sizeof(prime_field::field_element));
    __dst[1] = new prime_field::field_element[order];//(prime_field::field_element*)malloc(order * sizeof(prime_field::field_element));
    __dst[2] = new prime_field::field_element[order];
    twiddle_factor = new prime_field::field_element[order];
}

void delete_scratch_pad()
{
    delete[] __dst[0];
    delete[] __dst[1];
    delete[] twiddle_factor;
}

void fast_fourier_transform(const prime_field::field_element *coefficients, int coef_len, int order, prime_field::field_element root_of_unity, prime_field::field_element *result)
{
    prime_field::field_element rot_mul[62];
    //note: malloc and free will not call the constructor and destructor, not recommended unless for efficiency
    assert(sizeof(prime_field::field_element) * 2 == sizeof(__hhash_digest));
    //In sake of both memory and time efficiency, use the non-recursive version
    int lg_order = -1;
    rot_mul[0] = root_of_unity;
    for(int i = 0; i < 62; ++i)
    {
        if(i > 0)
            rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1];
        if((1LL << i) == order)
        {
            lg_order = i;
        }
    }
    int lg_coef = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == coef_len)
        {
            lg_coef = i;
        }
    }
    assert(lg_order != -1 && lg_coef != -1);
    // assert(rot_mul[lg_order].real == 1 && rot_mul[lg_order].img == 0);

    //we can merge both cases, but I just don't want to do so since it's easy to make mistake
    if(lg_coef > lg_order)
    {
        assert(false);
    }
    else
    {
        //initialize leaves
        int blk_sz = (order / coef_len);
        for(int j = 0; j < blk_sz; ++j)
        {
            for(int i = 0; i < coef_len; ++i)
            {
                __dst[lg_coef & 1][(j << lg_coef) | i] = coefficients[i];
            }
        }

        prime_field::field_element *x_arr = new prime_field::field_element[1 << lg_order];

        {
            //initialize leaves

            {
                for(int dep = lg_coef - 1; dep >= 0; --dep)
                {
                    int blk_size = 1 << (lg_order - dep);
                    int half_blk_size = blk_size >> 1;
                    int cur = dep & 1;
                    int pre = cur ^ 1;

                    prime_field::field_element x = prime_field::field_element(1);
                    x_arr[0] = prime_field::field_element(1);
                    for(int j = 1; j < blk_size; ++j)
                        x_arr[j] = x_arr[j - 1] * rot_mul[dep];
                    for(int k = 0; k < blk_size / 2; ++k)
                    {
                        int double_k = (k) & (half_blk_size - 1);
                        for(int j = 0; j < (1 << dep); ++j)
                        {
                            auto l_value = __dst[pre][double_k << (dep + 1) | j], r_value = x_arr[k] * __dst[pre][double_k << (dep + 1) | (1 << dep) | j];
                            __dst[cur][k << dep | j] = l_value + r_value;
                            __dst[cur][(k + blk_size / 2) << dep | j] = l_value - r_value;
                        }
                    }
                }
            }
        }
        delete[] x_arr;
    }

    for(int i = 0; i < order; ++i)
        result[i] = __dst[0][i];
}

void inverse_fast_fourier_transform(prime_field::field_element *evaluations, int coef_len, int order, prime_field::field_element root_of_unity, prime_field::field_element *dst)
{
    if(coef_len > order)
    {
        //more coefficient than evaluation
        fprintf(stderr, "Warning, Request do inverse fft with inefficient number of evaluations.");
        fprintf(stderr, "Will construct a polynomial with less order than required.");
        coef_len = order;
    }

    //assume coef_len <= order

    //subsample evalutions

    prime_field::field_element *sub_eval;
    bool need_free = false;
    if(coef_len != order)
    {
        need_free = true;
        sub_eval = (prime_field::field_element*)malloc(coef_len * sizeof(prime_field::field_element));
        for(int i = 0; i < coef_len; ++i)
        {
            sub_eval[i] = evaluations[i * (order / coef_len)];
        }
    }
    else
        sub_eval = evaluations;

    prime_field::field_element new_rou = prime_field::field_element(1);
    for(int i = 0; i < order / coef_len; ++i)
        new_rou = new_rou * root_of_unity;
    order = coef_len;

    prime_field::field_element inv_rou = prime_field::field_element(1), tmp = new_rou;
    int lg_order = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == order)
        {
            lg_order = i;
        }
    }
    int lg_coef = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == coef_len)
        {
            lg_coef = i;
        }
    }
    assert(lg_order != -1 && lg_coef != -1);

    for(int i = 0; i < lg_order; ++i)
    {
        inv_rou = inv_rou * tmp;
        tmp = tmp * tmp;
    }
    assert(inv_rou * new_rou == prime_field::field_element(1));

    fast_fourier_transform(sub_eval, order, coef_len, inv_rou, dst);

    if(need_free)
        free(sub_eval);

    prime_field::u256b two(2ull);
    prime_field::field_element inv_n = fast_pow(prime_field::field_element(order), prime_field::mod - two);
    assert(inv_n * prime_field::field_element(order) == prime_field::field_element(1));

    for(int i = 0; i < coef_len; ++i)
    {
        dst[i] = dst[i] * inv_n;
    }
}
