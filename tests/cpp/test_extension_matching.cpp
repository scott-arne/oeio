#include <gtest/gtest.h>

#include "oeio/oeio.h"

namespace oeio {
namespace test {

TEST(ExtensionMatchingTest, SimpleExtension) {
    const auto& registry = FormatRegistry::instance();
    const auto* handler = registry.lookup("test.sdf");
    ASSERT_NE(handler, nullptr);
    EXPECT_EQ(handler->info().name, "OEChem");
}

TEST(ExtensionMatchingTest, CompoundExtension) {
    const auto& registry = FormatRegistry::instance();
    const auto* handler = registry.lookup("test.sdf.gz");
    ASSERT_NE(handler, nullptr);
    EXPECT_EQ(handler->info().name, "OEChem");
}

TEST(ExtensionMatchingTest, PathWithDirectories) {
    const auto& registry = FormatRegistry::instance();
    const auto* handler = registry.lookup("/some/path/to/file.sdf");
    ASSERT_NE(handler, nullptr);
    EXPECT_EQ(handler->info().name, "OEChem");
}

TEST(ExtensionMatchingTest, LookupExtDirect) {
    const auto& registry = FormatRegistry::instance();
    const auto* handler = registry.lookup_ext(".sdf");
    EXPECT_NE(handler, nullptr);
}

TEST(ExtensionMatchingTest, LookupExtUnknown) {
    const auto& registry = FormatRegistry::instance();
    const auto* handler = registry.lookup_ext(".unknown");
    EXPECT_EQ(handler, nullptr);
}

TEST(ExtensionMatchingTest, OebGzExtension) {
    const auto& registry = FormatRegistry::instance();
    const auto* handler = registry.lookup("test.oeb.gz");
    ASSERT_NE(handler, nullptr);
    EXPECT_EQ(handler->info().name, "OEChem");
}

}  // namespace test
}  // namespace oeio
