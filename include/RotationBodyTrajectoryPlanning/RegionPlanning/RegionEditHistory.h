#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>

#include <cstddef>
#include <vector>

namespace smrobot::spray::rotationbody
{
    class RegionEditHistory
    {
    public:
        RegionEditHistory() = default;
        RegionEditHistory(SectionContour contour, RegionAssignment automaticAssignment);

        PlanningResult<void> resetAutomatic(
            SectionContour contour,
            RegionAssignment automaticAssignment);
        PlanningResult<void> applyRectangle(const YzRectangle& rectangle, RegionLabel label);

        bool canUndo() const noexcept;
        bool canRedo() const noexcept;
        bool undo() noexcept;
        bool redo() noexcept;
        void restoreAutomatic() noexcept;

        RegionAssignment resolved() const;
        const RegionAssignment& automatic() const noexcept;
        const std::vector<RegionOverrideCommand>& commands() const noexcept;
        std::size_t appliedCommandCount() const noexcept;

    private:
        SectionContour m_contour;
        RegionAssignment m_automatic;
        std::vector<RegionOverrideCommand> m_commands;
        std::size_t m_appliedCommandCount{ 0 };
    };
}
