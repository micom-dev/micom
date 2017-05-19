"""Test growth media agorithms."""

from fixtures import community
import numpy as np
import micom.media as media


class TestMinimalMedia:

    def test_linear(self, community):
        medium = media.minimal_medium(community, 0.8, 0.1)
        assert len(medium) <= 4
        assert all(medium > 1e-9)

    def test_mip(self, community):
        medium = media.minimal_medium(community, 0.8, 0.1,
                                      minimize_components=True)
        assert len(medium) <= 4
        assert all(medium > 1e-9)

        # Anaerobic growth
        medium = media.minimal_medium(community, 0.1, 0.1,
                                      minimize_components=True)
        assert len(medium) <= 3
        assert all(medium > 1e-9)

    def test_benchmark_linear(self, community, benchmark):
        benchmark(media.minimal_medium, community, 0.8, 0.1)

    def test_benchmark_mip(self, community, benchmark):
        benchmark(media.minimal_medium, community, 0.8, 0.1, True)
