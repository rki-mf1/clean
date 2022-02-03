include { download_host; check_own } from '../modules/get_host'

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow prepare_host {
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
    if ( params.own ) {
      check_own(ownFastaChannel)
      checkedOwn = check_own.out
    }
    else {
      checkedOwn = Channel.empty()
    }
  emit:
    host = host
    checkedOwn = checkedOwn
}