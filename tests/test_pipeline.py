# Copyright 2023, 2024 Julius Zimmermann
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for the volume conductor pipeline infrastructure."""

import numpy as np
import pytest

from ossdbs.utils import have_dielectric_properties_changed

pipeline_mod = pytest.importorskip(
    "ossdbs.fem.volume_conductor.pipeline",
    reason="pipeline module not yet implemented",
)
AnalysisPipeline = pipeline_mod.AnalysisPipeline
LazyFieldProvider = pipeline_mod.LazyFieldProvider
PipelineBuilder = pipeline_mod.PipelineBuilder
PipelineStage = pipeline_mod.PipelineStage
PreallocatedOutputArrays = pipeline_mod.PreallocatedOutputArrays

stages_mod = pytest.importorskip(
    "ossdbs.fem.volume_conductor.pipeline.stages",
    reason="pipeline.stages module not yet implemented",
)
ExportStage = stages_mod.ExportStage
FrequencySolverStage = stages_mod.FrequencySolverStage
MeshRefinementStage = stages_mod.MeshRefinementStage
PointEvaluationStage = stages_mod.PointEvaluationStage
SurfaceImpedanceSolverStage = stages_mod.SurfaceImpedanceSolverStage
TimeReconstructionStage = stages_mod.TimeReconstructionStage


class TestPipelineStage:
    """Tests for PipelineStage ABC."""

    def test_abstract_methods(self):
        """Test that PipelineStage cannot be instantiated directly."""
        with pytest.raises(TypeError):
            PipelineStage()

    def test_concrete_stage_has_name(self):
        """Test that concrete stages have a name property."""
        stage = MeshRefinementStage()
        assert stage.name == "MeshRefinement"

        stage = FrequencySolverStage()
        assert stage.name == "FrequencySolver"

        stage = PointEvaluationStage()
        assert stage.name == "PointEvaluation"

        stage = TimeReconstructionStage()
        assert stage.name == "TimeReconstruction"

        stage = ExportStage()
        assert stage.name == "Export"


class TestAnalysisPipeline:
    """Tests for AnalysisPipeline."""

    def test_empty_pipeline(self):
        """Test creating an empty pipeline."""
        pipeline = AnalysisPipeline()
        assert len(pipeline.stages) == 0

    def test_add_stage(self):
        """Test adding stages to pipeline."""
        pipeline = AnalysisPipeline()
        stage = MeshRefinementStage()
        pipeline.add_stage(stage)
        assert len(pipeline.stages) == 1
        assert pipeline.stages[0].name == "MeshRefinement"

    def test_add_stage_chaining(self):
        """Test that add_stage returns self for chaining."""
        pipeline = AnalysisPipeline()
        result = pipeline.add_stage(MeshRefinementStage())
        assert result is pipeline

    def test_insert_stage(self):
        """Test inserting a stage at a specific position."""
        pipeline = AnalysisPipeline()
        pipeline.add_stage(MeshRefinementStage())
        pipeline.add_stage(ExportStage())
        pipeline.insert_stage(1, FrequencySolverStage())

        assert len(pipeline.stages) == 3
        assert pipeline.stages[0].name == "MeshRefinement"
        assert pipeline.stages[1].name == "FrequencySolver"
        assert pipeline.stages[2].name == "Export"

    def test_remove_stage(self):
        """Test removing a stage by name."""
        pipeline = AnalysisPipeline()
        pipeline.add_stage(MeshRefinementStage())
        pipeline.add_stage(FrequencySolverStage())
        pipeline.add_stage(ExportStage())

        pipeline.remove_stage("FrequencySolver")

        assert len(pipeline.stages) == 2
        assert pipeline.stages[0].name == "MeshRefinement"
        assert pipeline.stages[1].name == "Export"

    def test_stages_returns_copy(self):
        """Test that stages property returns a copy."""
        pipeline = AnalysisPipeline()
        pipeline.add_stage(MeshRefinementStage())

        stages = pipeline.stages
        stages.clear()

        assert len(pipeline.stages) == 1


class TestPipelineBuilder:
    """Tests for PipelineBuilder."""

    def test_default_pipeline(self):
        """Test creating a default pipeline."""
        pipeline = PipelineBuilder.default_pipeline()
        assert len(pipeline.stages) == 5

    def test_builder_chaining(self):
        """Test builder method chaining."""
        pipeline = (
            PipelineBuilder()
            .with_mesh_refinement()
            .with_frequency_solver()
            .with_export()
            .build()
        )
        assert len(pipeline.stages) == 3


class TestLazyFieldProvider:
    """Tests for LazyFieldProvider optimization."""

    def test_field_not_computed_on_init(self):
        """Test that field is not computed during initialization."""
        # We can't test with real ngsolve objects without a mesh,
        # but we can test the lazy behavior
        provider = LazyFieldProvider(None)
        assert provider._field_computed is False
        assert provider._field is None

    def test_invalidate(self):
        """Test invalidating cached field."""
        provider = LazyFieldProvider(None)
        provider._field_computed = True
        provider._field = "cached"

        provider.invalidate()

        assert provider._field_computed is False
        assert provider._field is None

    def test_potential_setter_invalidates(self):
        """Test that setting potential invalidates cached field."""
        provider = LazyFieldProvider(None)
        provider._field_computed = True
        provider._field = "cached"

        provider.potential = "new_potential"

        assert provider._field_computed is False
        assert provider._field is None
        assert provider.potential == "new_potential"


class TestPreallocatedOutputArrays:
    """Tests for PreallocatedOutputArrays optimization."""

    def test_array_shapes(self):
        """Test that arrays have correct shapes."""
        n_points = 1000
        arrays = PreallocatedOutputArrays(n_points, complex)

        assert arrays.potential.shape == (n_points,)
        assert arrays.field.shape == (n_points, 3)
        assert arrays.potential.dtype == complex
        assert arrays.field.dtype == complex

    def test_float_dtype(self):
        """Test with float dtype."""
        n_points = 500
        arrays = PreallocatedOutputArrays(n_points, float)

        assert arrays.potential.dtype == float
        assert arrays.field.dtype == float

    def test_reset(self):
        """Test resetting arrays to zero."""
        arrays = PreallocatedOutputArrays(10, float)
        arrays.potential[:] = 1.0
        arrays.field[:] = 1.0

        arrays.reset()

        np.testing.assert_array_equal(arrays.potential, 0.0)
        np.testing.assert_array_equal(arrays.field, 0.0)


class TestMeshRefinementStage:
    """Tests for MeshRefinementStage."""

    def test_stage_name(self):
        """Test stage name."""
        stage = MeshRefinementStage()
        assert stage.name == "MeshRefinement"

    def test_init_with_parameters(self):
        """Test initialization with parameters."""
        stage = MeshRefinementStage(
            material_steps=2,
            adaptive_settings={
                "Active": True,
                "ErrorTolerance": 1.0,
                "MaxIterations": 5,
            },
        )
        assert stage._material_steps == 2
        assert stage._adaptive_settings["Active"] is True

    def test_validate_adaptive_settings_missing_keys(self):
        """Test validation fails with missing keys."""
        stage = MeshRefinementStage()
        with pytest.raises(ValueError, match="ErrorTolerance"):
            stage._validate_adaptive_settings({"Active": True})


class TestTimeReconstructionStage:
    """Tests for TimeReconstructionStage."""

    def test_stage_name(self):
        """Test stage name."""
        stage = TimeReconstructionStage()
        assert stage.name == "TimeReconstruction"


class TestPointEvaluationStage:
    """Tests for PointEvaluationStage."""

    def test_stage_name(self):
        """Test stage name."""
        stage = PointEvaluationStage()
        assert stage.name == "PointEvaluation"


class TestSurfaceImpedanceSolverStage:
    """Tests for SurfaceImpedanceSolverStage."""

    def test_stage_name(self):
        """Test stage name."""
        stage = SurfaceImpedanceSolverStage()
        assert stage.name == "SurfaceImpedanceSolver"

    def test_default_parameters(self):
        """Test default parameter values."""
        stage = SurfaceImpedanceSolverStage()
        assert stage._estimate_currents is False
        assert stage._dielectric_threshold == 0.01

    def test_custom_parameters(self):
        """Test custom parameter values."""
        stage = SurfaceImpedanceSolverStage(
            estimate_currents=True,
            dielectric_threshold=0.05,
        )
        assert stage._estimate_currents is True
        assert stage._dielectric_threshold == 0.05


class TestHaveDielectricPropertiesChanged:
    """Tests for have_dielectric_properties_changed utility."""

    def test_returns_true_for_significant_change(self):
        """Test returns True when change exceeds threshold."""

        class MockModel:
            def conductivity(self, omega):
                # Returns different values at different frequencies
                if omega < 1000:
                    return 1.0
                return 2.0

        props = {"material": MockModel()}
        result = have_dielectric_properties_changed(
            props, is_complex=False, old_freq=100, new_freq=200, threshold=0.01
        )
        assert result is True

    def test_returns_false_for_small_change(self):
        """Test returns False when change is below threshold."""

        class MockModel:
            def conductivity(self, omega):
                # Returns nearly constant value
                return 1.0 + omega * 1e-10

        props = {"material": MockModel()}
        result = have_dielectric_properties_changed(
            props, is_complex=False, old_freq=100, new_freq=200, threshold=0.01
        )
        assert result is False

    def test_complex_conductivity(self):
        """Test with complex conductivity values."""

        class MockModel:
            def conductivity(self, omega):
                return 1.0 + 0.1j * (omega / 1000)

        props = {"material": MockModel()}
        result = have_dielectric_properties_changed(
            props, is_complex=True, old_freq=100, new_freq=1000, threshold=0.01
        )
        # The imaginary part changes significantly
        assert result is True

    def test_empty_properties(self):
        """Test with empty dielectric properties."""
        result = have_dielectric_properties_changed(
            {}, is_complex=False, old_freq=100, new_freq=200, threshold=0.01
        )
        # Should return False since no materials changed
        assert result is False
