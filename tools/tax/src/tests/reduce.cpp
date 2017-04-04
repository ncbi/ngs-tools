#include "tests.h"

TEST(reduce_1) {
    vector<int> v = {1, 2};
    v = reduce::median_filter(v, 2);
    equal(v.size(), 2);
    equal(v[0], 1);
    equal(v[1], 2);
};

TEST(reduce_2) 	{
    vector<int> v = {1, 2, 3, 4};
    v = reduce::median_filter(v, 2);
    equal(v.size(), 4);
    equal(v[0], 1);
    equal(v[1], 2);
    equal(v[2], 3);
    equal(v[3], 4);
};

TEST(reduce_3) {
    vector<int> v = {1, 2, 3, 4, 5};
    v = reduce::median_filter(v, 2);
    equal(v.size(), 5);
    equal(v[0], 1);
    equal(v[1], 2);
    equal(v[2], 3);
    equal(v[3], 4);
    equal(v[4], 5);
};

TEST(reduce_4) {													 //10			// 14
    vector<int> _v = {1, 2, 3, -40, 5, 6, 7, 8, 9, 10, 40, 11, 244, 5, 2, -10, 4, 5};
    auto v = reduce::median_filter(_v, 2);
    equal(v.size(), 18);
    equal(v[0], 1);
    equal(v[1], 2);
    equal(v[2], 2);
    equal(v[3], 3);
    equal(v[4], 5);
    equal(v[5], 6);
    equal(v[6], 7);
    equal(v[7], 8);
    equal(v[8], 9);
    equal(v[9], 10);
    equal(v[10], 11); // 40
    equal(v[11], 11);
    equal(v[12], 11);
    equal(v[13], 5);

    equal(v[14], 4);
    equal(v[15], 4);
    equal(v[16], 4);
    equal(v[17], 5);
};

TEST_MAIN();
