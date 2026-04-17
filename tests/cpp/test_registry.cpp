#include <gtest/gtest.h>

#include "oeio/oeio.h"

// Force-link the OEChem handler when using static libraries
namespace oeio {
extern void oeio_force_link_oechem_handler();
}

namespace {
struct ForceLink {
    ForceLink() { oeio::oeio_force_link_oechem_handler(); }
} force_link_instance;
}

namespace oeio {
namespace test {

TEST(FormatRegistryTest, SingletonReturnsConsistentInstance) {
    const auto& instance1 = FormatRegistry::instance();
    const auto& instance2 = FormatRegistry::instance();
    EXPECT_EQ(&instance1, &instance2);
}

TEST(FormatRegistryTest, OEChemHandlerRegisteredByDefault) {
    const auto& registry = FormatRegistry::instance();
    const auto* handler = registry.lookup("test.sdf");
    EXPECT_NE(handler, nullptr);
}

TEST(FormatRegistryTest, LookupReturnsNullForUnknownExtension) {
    const auto& registry = FormatRegistry::instance();
    const auto* handler = registry.lookup("test.xyz123");
    EXPECT_EQ(handler, nullptr);
}

TEST(FormatRegistryTest, FormatsReturnsNonEmpty) {
    const auto& registry = FormatRegistry::instance();
    auto formats = registry.formats();
    EXPECT_GT(formats.size(), 0);
}

TEST(FormatRegistryTest, FormatsContainsOEChem) {
    const auto& registry = FormatRegistry::instance();
    auto formats = registry.formats();
    bool found_oechem = false;
    for (const auto& fmt : formats) {
        if (fmt.name == "OEChem") {
            found_oechem = true;
            break;
        }
    }
    EXPECT_TRUE(found_oechem);
}

}  // namespace test
}  // namespace oeio
