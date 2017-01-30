#ifndef FEATURES_H_INCLUDED
#define FEATURES_H_INCLUDED

#include <map>
#include <cmath>
#include <iostream>
#include <iomanip>

typedef int LowerLevelFeatureIndex; // actually char for first iteration
typedef double LowerLevelFeatureCount;

struct Features : public std::map<LowerLevelFeatureIndex, LowerLevelFeatureCount>
{
    double norm;
    Features() : norm(0){}
};

struct FeatureOperations
{
    static double norm(const Features &fs)
    {
        double s = 0;
        for (auto &f : fs)
            s += f.second * f.second;

        return std::sqrt(s);
    }

    static void add(Features &target, const Features &fs)
    {
        for (auto &f : fs)
            target[f.first] += f.second;

        target.norm = norm(target);
    }

    static double mul(const Features &a, const Features &b)
    {
        double s = 0;

        for (auto &f : a)
        {
            auto x = f.first;
            auto it = b.find(x);
            if (it != b.end())
                s += f.second * it->second;
        }

//        std::cout << "mul: " << s << std::endl;

        return s;
    }

    static double get_similarity(const Features &a, const Features &b)
    {
        auto norm_mul = a.norm * b.norm;
        if (norm_mul <= 0.00001) // todo: think
            return 0;
            //throw std::runtime_error("get_similarity norm_mul <= 0");

  //      std::cout << "norm mul: " << norm_mul << std::endl;

        double sim = mul(a, b)/norm_mul;
        if (sim > 1.0000001)
            throw std::runtime_error("sim > 1.0000001");

        return std::min(1.0, sim);
    }

    static Features calculate_features(const char *s, int kmer_len)
    {   
        Features fs;
        for (int i = 0; i < kmer_len; i++)
            fs[s[i]] ++;

        fs.norm = norm(fs);
        return fs;
    }
};

template <class IDCast = LowerLevelFeatureIndex>
static void print_features(const Features &fs)
{
    std::cout << "norm " << fs.norm << " ";
    for (auto &f : fs)
        std::cout << IDCast(f.first) << " : " << size_t(f.second + 0.5) << '\t';
//        std::cout << IDCast(f.first) << " : " << std::setprecision(2) << f.second << '\t';

    std::cout << std::endl;
}

#endif