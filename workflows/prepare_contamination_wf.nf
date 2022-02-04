include { download_host; check_own; concat_contamination } from '../modules/prepare_contamination'

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow prepare_contamination {
  take:
    nanoControlFastaChannel
    illuminaControlFastaChannel
    rRNAChannel

  main:
    prepare_auto_host()
    prepare_own_host()

    contamination_collection = prepare_auto_host.out.collect()
          .mix(nanoControlFastaChannel)
          .mix(illuminaControlFastaChannel)
          .mix(prepare_own_host.out)
          .mix(rRNAChannel).collect()
    concat_contamination(params.name, contamination_collection)

  emit:
    concat_contamination.out.fa
}

workflow prepare_auto_host {
  main:
    if ( params.host ) {
      if ( params.cloudProcess ) {
        host_preload = file("${params.databases}/hosts/${params.host}.fa.gz")
        if ( host_preload.exists() ) {
          host = Channel.fromPath(host_preload)
        } else {
          download_host(hostNameChannel)
          host = download_host.out
        }
      } else {
        download_host(hostNameChannel)
        host = download_host.out
      }
    }
    else {
      host = Channel.empty()
    }
  emit:
    host
}

workflow prepare_own_host {
  main:
    if ( params.own ) {
      check_own(ownFastaChannel)
      checkedOwn = check_own.out
    }
    else {
      checkedOwn = Channel.empty()
    }
  emit:
    checkedOwn
}