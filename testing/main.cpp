#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include "test_utils.hpp"
int main( int argc, char **argv)
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest( TestUtils::suite() );
  runner.run();
  return 0;
}