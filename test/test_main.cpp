/**
 * Test Main Runner
 * Runs all registered test cases
 */

#include "test_framework.h"

int main() {
    return TestRunner::instance().run_all();
}
