# pylint: disable=R0913

""" alignment_runner """

import subprocess


class AlignmentRunner:
    def __init__(self, cpus, ref_index, sample_id="sample_x", blocksize=10_000_000):
        self.sample_id = sample_id
        self.aligner = None
        self.aligner_args = [
            f"-t {cpus}",
            f"-K {blocksize}",
            ref_index,
        ]

    def generate_aligner_call(self, input_files=None, single_end_reads=False, alignment_file=None):
        read_group = f"'@RG\\tID:{1 if single_end_reads else 2}\\tSM:{self.sample_id}'"
        aligner_call = None
        if input_files is not None:
            aligner_call = f"{self.aligner} -R {read_group} {' '.join(self.aligner_args)} {' '.join(input_files)}"
            if alignment_file is not None:
                aligner_call = " | ".join([
                    aligner_call,
                    f"tee -a {alignment_file}"
                ])
        return aligner_call

    def run(self, input_files, single_end_reads=False, alignment_file=None):
        aligner_call = self.generate_aligner_call(
            input_files=input_files,
            single_end_reads=single_end_reads,
            alignment_file=alignment_file
        )

        # pylint: disable=R1732
        read_processing_proc = subprocess.Popen(aligner_call, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        return read_processing_proc, aligner_call  # .stdout

        # with subprocess.Popen(aligner_call, shell=True, stdout=subprocess.PIPE) as read_processing_proc:
        #    return read_processing_proc.stdout

        # try:
        #     with subprocess.Popen(
        #         commands, shell=True, stdout=subprocess.PIPE
        #     ) as read_processing_proc:
        #         profiler.count_alignments(
        #             read_processing_proc.stdout,
        #             aln_format="sam",
        #             min_identity=min_identity,
        #             min_seqlen=min_seqlen,
        #         )
        # except Exception as err:
        #     if isinstance(err, ValueError) and str(err).strip() == "file does not contain alignment data":
        #         logger.error(f"Failed to align. Is `{self.aligner}` installed and on the path?")
        #         sys.exit(1)

        #     logger.error("Caught some exception:")
        #     logger.error("%s", err)
        #     raise Exception from err


class BwaMemRunner(AlignmentRunner):
    def __init__(self, cpus, ref_index, sample_id="sample_x", blocksize=10_000_000):
        super().__init__(cpus, ref_index, sample_id=sample_id, blocksize=blocksize)
        self.aligner = "bwa mem"
        self.aligner_args += [
            "-v 1",
            "-a",
        ]


class Minimap2Runner(AlignmentRunner):
    def __init__(self, cpus, ref_index, sample_id="sample_x", blocksize=10_000_000):
        super().__init__(cpus, ref_index, sample_id=sample_id, blocksize=blocksize)
        self.aligner = "minimap2"
        self.aligner_args += [
            "--sam-hit-only",
            "-x sr",
            "--secondary=yes",
            "-a",
            f"--split-prefix {sample_id}_split",
        ]
