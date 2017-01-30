#include "feature_tree_builder.h"
#include "feature_tree_io.h"

using namespace std;

template <class A, class B>
void equal(A a, B b)
{
	if (a == b)
		cout << "ok" << endl;
	else
	{
		cout << "FAIL: " << a << " != " << b << endl;
//		exit(1);
	}
}

Features fs(const string &s)
{
    return FeatureOperations::calculate_features(&s[0], s.length());
}

void find_closest(const char *nine_mer, const Tree &tree_3mers, const Tree &tree_9mers)
{
    SimilarityOf3Mers sim3(tree_3mers);

    auto fs = calculate_features_for_3mers(nine_mer, 9, tree_3mers);
//    print_features(fs);

    auto closest = find_closest_final_node(tree_9mers, fs, sim3);
    cout << "closest id: " << closest.id << " similarity: " << closest.similarity << " feature norm " << tree_9mers.get_node(closest.id).features.norm << endl;
}

void find_closest(const char *kmer, const Tree &tree_3mers, const Tree &tree_9mers, const Tree &tree_27mers)
{
    auto fs = calculate_features_for_9mers(kmer, 27, tree_3mers, tree_9mers);
//    print_features(fs);

    SimilarityOfNMers sim9(tree_9mers);
    auto closest = find_closest_final_node(tree_27mers, fs, sim9);
    if (!tree_27mers.valid_node(closest.id))
    {
        cout << "not found" << endl;
        return;
    }

    cout << "closest id: " << closest.id << " similarity: " << closest.similarity << " feature norm " << tree_27mers.get_node(closest.id).features.norm << endl;
    tree_27mers.for_node_and_parents_do(closest.id, [&](const Tree::Node &node) 
    {
        cout << " " << node.id;
    });

    cout << endl;
}

int main(int argc, char const *argv[])
{
    try
    {
        {
            Tree tree_3mers; 
            TreeIO::load_tree(tree_3mers, "./tests/tree.3mer.0.5.save");

            Tree tree_9mers; 
            TreeIO::load_tree(tree_9mers, "./tests/tree.3mer_x_3.0.12.9.7.save");

            Tree tree_27mers; 
            TreeIO::load_tree(tree_27mers, "./tests/tree_27mer.0.17.27.25.overnight.12.save");

            const char *EBOLA = "TGAAATTGTTACTGTAATCATACCTGGTTTGTTTCAGAGCCATATCACCAAGATAGAGAACAACCTAGGTCTCCGGAGGGGGCAAGGGCATCAGTGTGCTCAGTTGAAAATCCCTTGTCAACATCTAGGCCTTATCACATCACAAGTTCCGCCTTAAACTCTGCAGGGTGATCC";
            for (int i=0; i < 30; i++)
            {
                if (i % 3 == 0)
                    cout << endl;
                cout << "i " << i << " ";
                find_closest(&EBOLA[i], tree_3mers, tree_9mers, tree_27mers);
            }

            return 0;
        }


#if 0
        {
            Tree tree_3mers; 
            TreeIO::load_tree(tree_3mers, "./tests/tree.3mer.0.5.save");

            Tree tree_9mers; 
            TreeIO::load_tree(tree_9mers, "./tests/tree.3mer_x_3.0.12.9.7.save");

            auto fa = calculate_features_for_9mers("GAGCCATATCACCAAGATAGAGAACAA", 21, tree_3mers, tree_9mers);
            auto fb = calculate_features_for_9mers("AATAGCAACATTATTGTTAAAGGACAG", 21, tree_3mers, tree_9mers);
            {
                SimilarityOf1Mers sim1;
                cout << "sim 1: " << sim1(fa, fb) << endl;
            }

            {
                SimilarityOfNMers simn(tree_9mers);
                cout << "sim n: " << simn(fa, fb) << endl;
            }
#if 0

            const char *EBOLA = "TGAAATTGTTACTGTAATCATACCTGGTTTGTTTCAGAGCCATATCACCAAGATAGAGAACAACCTAGGTCTCCGGAGGG";
            for (int i=0; i < 30; i++)
            {
                if (i % 3 == 0)
                    cout << endl;
                cout << "i " << i << " ";
                find_closest(&EBOLA[i], tree_3mers, tree_9mers);
            }
#endif
//            find_closest("TGAAATTGT", tree_3mers, tree_9mers);
  //          find_closest("GAAATTGTT", tree_3mers, tree_9mers);
    //        find_closest("AAATTGTTA", tree_3mers, tree_9mers);
            return 0;
        }
#endif
        {
            Tree tree_3mers; 
            TreeIO::load_tree(tree_3mers, "./tests/tree.3mer.0.5.save");
#if 0
            auto fs = calculate_features_for_3mers("ACACTGAAG", 9, tree_3mers);
            cout << "features size: " << fs.size() << endl;
            cout << "norm " << fs.norm << " ";
            print_features(fs);
            cout << "-------------" << endl;
            cout << "extended:" << endl;
            auto extended = extend_features_to_parents(fs, tree_3mers);
            cout << "norm " << extended.norm << " ";
            print_features(extended);

            {
                auto fa = calculate_features_for_3mers("ACACTGAAG", 9, tree_3mers);
                auto fb = calculate_features_for_3mers("CGCTAAAGA", 9, tree_3mers);
                {
                    SimilarityOf1Mers sim1;
                    cout << "sim 1: " << sim1(fa, fb) << endl;
                }

                {
                    SimilarityOf3Mers sim3(tree_3mers);
                    cout << "sim 3: " << sim3(fa, fb) << endl;
                }
            }

            {
                auto fa = calculate_features_for_3mers("ACACTGAAG", 9, tree_3mers);
                auto fb = calculate_features_for_3mers("ACACTGAAG", 9, tree_3mers);
                {
                    SimilarityOf1Mers sim1;
                    cout << "sim 1: " << sim1(fa, fb) << endl;
                }

                {
                    SimilarityOf3Mers sim3(tree_3mers);
                    cout << "sim 3: " << sim3(fa, fb) << endl;
                }
            }

            {
                auto fa = calculate_features_for_3mers("ACACTGAGG", 9, tree_3mers);
                cout << "fa ";
                print_features(fa);
                auto levels = SimilarityOf3Mers::get_levels(fa, tree_3mers);
                for (auto &l : levels)
                {
                    cout << l.first << " ";
                    print_features(l.second);
                }
            }

            {
                auto fa = calculate_features_for_3mers("ACA", 3, tree_3mers);
                cout << "fa ";
                print_features(fa);
                auto levels = SimilarityOf3Mers::get_levels(fa, tree_3mers);
                for (auto &l : levels)
                {
                    cout << l.first << " ";
                    print_features(l.second);
                }
            }
#endif

            {
                cout << "height of 0 " << tree_3mers.get_height_level(0) << endl;
                cout << "height of 7 " << tree_3mers.get_height_level(7) << endl;

                auto fa = calculate_features_for_3mers("TGAAATTGT", 9, tree_3mers);
                auto fb = calculate_features_for_3mers("AACCTAGGT", 9, tree_3mers);
                cout << "fa " << fa.norm << " ";
                print_features(fa);
                cout << "fb " << fb.norm << " ";
                print_features(fb);
                {
                    SimilarityOf1Mers sim1;
                    cout << "sim 1: " << sim1(fa, fb) << endl;
                }

                {
                    SimilarityOf3Mers sim3(tree_3mers);
                    cout << "sim 3: " << sim3(fa, fb) << endl;
                }

                //cout << "extended: fa_ex "; 
                //auto fa_ex = extend_features_to_parents(fa, tree_3mers);
                //print_features(fa_ex);

                //cout << "extended: fb_ex ";
                //auto fb_ex = extend_features_to_parents(fb, tree_3mers);
                //print_features(fb_ex);
            }

        }

        return 0;

        {
            SimilarityOf1Mers sim;
            TreeBuilder<SimilarityOf1Mers> b(sim);

            b.insert(fs("good"));
            TreeIO::print_tree<char>(b.tree, b.tree.head_id);
            cout << "-------------" << endl;

            b.insert(FeatureOperations::calculate_features("bad", 3));
            TreeIO::print_tree<char>(b.tree, b.tree.head_id);
            cout << "-------------" << endl;

            b.insert(FeatureOperations::calculate_features("bod", 3));
            TreeIO::print_tree<char>(b.tree, b.tree.head_id);
            cout << "-------------" << endl;

            b.insert(FeatureOperations::calculate_features("dad", 3));
            b.insert(FeatureOperations::calculate_features("bid", 3));
            TreeIO::print_tree<char>(b.tree, b.tree.head_id);
            cout << "-------------" << endl;


            {
                auto c = find_closest_final_node(b.tree, fs("good"), sim);
                cout << c.id << " " << c.similarity << endl;
            }
            {
                auto c = find_closest_final_node(b.tree, fs("bild"), sim);
                cout << c.id << " " << c.similarity << endl;
            }
            {
                auto c = find_closest_final_node(b.tree, fs("bad"), sim);
                cout << c.id << " " << c.similarity << endl;
            }
        }

        return 0;
    }
    catch ( exception & x )
    {
        cerr << x.what() << endl;
		cerr << "exit 3" << endl;
		return 3;
    }
    catch ( string & x )
    {
        cerr << x << endl;
		cerr << "exit 4" << endl;
		return 4;
    }
    catch ( const char * x )
    {
        cerr << x << endl;
		cerr << "exit 5" << endl;
		return 5;
    }
    catch ( ... )
    {
        cerr << "unknown exception" << endl;
//		cerr << "exit 6" << endl;
		return 6;
    }
}
