/**
 * Minimal Test Framework Header
 * Simple assertion-based testing without external dependencies
 */

#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

#include <iostream>
#include <string>
#include <functional>
#include <vector>
#include <iomanip>
#include <sstream>

class TestRunner {
public:
    struct TestCase {
        std::string name;
        std::string section;
        std::function<void()> test_func;
    };

    static TestRunner& instance() {
        static TestRunner runner;
        return runner;
    }

    static void add_test(const std::string& name, const std::string& section, std::function<void()> func) {
        instance().tests.push_back({name, section, func});
    }

    int run_all() {
        int passed = 0;
        int failed = 0;
        int test_num = 0;
        
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "Running " << tests.size() << " test cases...\n";
        std::cout << std::string(70, '=') << "\n\n";

        for (const auto& test : tests) {
            test_num++;
            try {
                std::cout << "[" << std::setw(3) << test_num << "/" << std::setw(3) << tests.size() << "] ";
                std::cout << test.name << " :: " << test.section << " ... ";
                std::cout.flush();
                test.test_func();
                std::cout << "\033[32mPASS\033[0m\n";
                passed++;
            } catch (const std::exception& e) {
                std::cout << "\033[31mFAIL\033[0m\n";
                std::cout << "         Error: " << e.what() << "\n";
                failed++;
            }
        }

        std::cout << "\n" << std::string(70, '-') << "\n";
        std::cout << "Results: " << passed << " passed, " << failed << " failed\n";
        std::cout << std::string(70, '=') << "\n\n";

        return failed == 0 ? 0 : 1;
    }

private:
    std::vector<TestCase> tests;
};

#define REQUIRE(condition) \
    do { \
        if (!(condition)) { \
            std::ostringstream oss; \
            oss << "Assertion failed at " << __FILE__ << ":" << __LINE__ << "\n"; \
            oss << "  Condition: " << #condition; \
            throw std::runtime_error(oss.str()); \
        } \
    } while(0)

#define REQUIRE_EQ(a, b) \
    do { \
        if ((a) != (b)) { \
            std::ostringstream oss; \
            oss << "Assertion failed at " << __FILE__ << ":" << __LINE__ << "\n"; \
            oss << "  Expected: " << (b) << "\n"; \
            oss << "  Got: " << (a); \
            throw std::runtime_error(oss.str()); \
        } \
    } while(0)

#define PP_CAT(a, b) PP_CAT_I(a, b)
#define PP_CAT_I(a, b) a ## b

#define TEST_CASE(name, tags) \
    static void PP_CAT(test_impl_, __LINE__)(); \
    namespace { \
        struct PP_CAT(TestRegister_, __LINE__) { \
            PP_CAT(TestRegister_, __LINE__)() { \
                TestRunner::add_test(name, tags, &PP_CAT(test_impl_, __LINE__)); \
            } \
        } PP_CAT(register_instance_, __LINE__); \
    } \
    static void PP_CAT(test_impl_, __LINE__)()

#endif // TEST_FRAMEWORK_H
